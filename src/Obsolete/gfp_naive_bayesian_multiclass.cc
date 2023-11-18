#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <map>
#include <vector>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Utilities/GFP_Tools/gfp.h"

using std::map;
using std::pair;
using std::vector;

using std::cerr;
using std::endl;

const char * prog_name = nullptr;


static int strip_leading_zeros = 0;
static int ignore_missing_activity = 0;
static int activity_column = 1;   // by default identifier in column 0
typedef vector<int> Class_IDs;
typedef map<IWString,Class_IDs> Activity;
static Activity activity;
static IW_STL_Hash_Map_int class_name_to_id;
static vector<IWString> class_names;
static int number_classes = 0;

static int bit_count_cutoff = 1;
static int include_lower_counts = 1;

typedef pair<unsigned int, unsigned int> Bit_And_Count;
typedef vector<double> Frequencies;
typedef pair<int,Bit_And_Count> FP_And_Bit_And_Count;
typedef map<FP_And_Bit_And_Count,Frequencies> Bits_Profiler;
static Bits_Profiler bits_profiler;

static int number_fixed=0;
static int number_sparse=0;

static Frequencies training_set;

#define HEADER_RECORD "# written by gfp_naive_bayesian_multiclass v2"
#define COUNT_FIXED "count_fixed"
#define COUNT_SPARSE "count_sparse"
#define CLASS_NAMES  "class_labels"
#define BIT_COUNT_CUTOFF "bit_count_cutoff"
#define INCLUDE_LOWER_COUNTS "include_lower_counts"

static int verbose = 0;

static int train_mode=1;
static int header_in_activity_file = 0;
static int min_support_level = 0;
static int write_probabilities_in_01_range = 0;
static double laplace_correction_factor=1.0;

static int
read_activity_data_record (const const_IWSubstring & buffer,
                           Activity & activity)
{
  IWString id, token;
  int i = 0;

  char separator = ' ';

  if (buffer.nextword(id, i) && buffer.nextword(token, i))
    ;
  else if (buffer.nextword(id, i, '\t') && buffer.nextword(token, i, '\t'))
    separator = '\t';
  else if (buffer.nextword(id, i, ',') && buffer.nextword(token, i, ','))
    separator = ',';
  else
  {
    cerr << "Cannot separate into identifier and activity '" << buffer << "'\n";
    return 0;
  }

  if (1 != activity_column)
  {
    if (! buffer.word(activity_column, token, separator))
    {
      cerr << "Cannot extract column '" << (activity_column + 1) << " from record\n";
      return 0;
    }
  }

  if (strip_leading_zeros)
    id.remove_leading_chars('0');

  if (!class_name_to_id.contains(token))
  {
    class_name_to_id[token] = class_names.size();
    class_names.push_back(token);
  }

  if (activity.find(id) != activity.end())
  {
    activity[id].push_back(class_name_to_id[token]);
  }
  else
  {
    activity[id]=Class_IDs();
    activity[id].push_back(class_name_to_id[token]);
  }

//cerr << "for id '" << id << "' value '" << activity[id] << "', token '" << token << "'\n";

  return 1;
}

static int
read_activity_data (iwstring_data_source & input,
                    Activity & activity)
{
  const_IWSubstring buffer;
  while (input.next_record(buffer))
  {
    if (buffer.starts_with('#'))
      continue;
    if(header_in_activity_file && 1== input.lines_read()) //header line
      continue;
    if (! read_activity_data_record(buffer, activity))
    {
      cerr << "Cannot read activity data, line " << input.lines_read() << endl;
      cerr << "'" << buffer << "'\n";
      return 0;
    }
  }

  number_classes = class_names.size();
  training_set.resize(number_classes);
  return activity.size();
}

static int
read_activity_data (const const_IWSubstring & fname,
                    Activity & activity)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(1);

  return read_activity_data (input, activity);
}

static int
find_identifier_in_activity_hash (const Activity & activity,
                                  const IWString & id,
                                  Class_IDs & activity_for_id)
{
  Activity::const_iterator f = activity.find(id);

  if (f != activity.end())
  {
    activity_for_id = (*f).second;

    return 1;
  }

  if (! strip_leading_zeros)
    return 0;

  if (! id.starts_with('0'))
    return 0;

  IWString tmp(id);

  tmp.remove_leading_chars('0');

  f = activity.find(tmp);

  if (f == activity.end())
    return 0;

  activity_for_id = (*f).second;

  return 1;
}

static int
read_bits_profiler_record (const const_IWSubstring & buffer,
                           Bits_Profiler & bits_profiler,
                           const int & fpi)
{
  int i = 0;
  const_IWSubstring token;

  if (! buffer.nextword(token, i))    // skip blank lines
    return 0;

  unsigned int b;
  if (! token.numeric_value(b))
  {
    cerr << "Invalid bit number '" << buffer << "'\n";
    return 0;
  }
  
  if (! buffer.nextword(token, i))    // skip blank lines
    return 0;

  unsigned int c;
  if (! token.numeric_value(c))
  {
    cerr << "Invalid bit count '" << buffer << "'\n";
    return 0;
  }
 
  Frequencies freq(number_classes);
  for (int nc =0; nc < number_classes; ++nc)
  {
    if (! buffer.nextword(token, i))
    {
      cerr << "bits profiler records must contain at least " << number_classes<< " tokens '" << buffer << "'\n";
      return 0;
    }

    double num;

    if (! token.numeric_value(num) || num < 0.0)
    {
      cerr << "Frequencies must be positive real numbers '" << buffer << "'\n";
      return 0;
    }
    freq[nc] = num;
  }

  FP_And_Bit_And_Count fabc(fpi,Bit_And_Count(b,c));
  bits_profiler[fabc] = freq;

  return 1;
}

static int
read_bits_profiler(Bits_Profiler & bits_profiler,iwstring_data_source & input)
{
  const_IWSubstring buffer;

  if (! input.next_record(buffer))
  {
    cerr << "read_bits_profiler:cannot read header\n";
    return 0;
  }

  buffer.strip_trailing_blanks();

  if (HEADER_RECORD != buffer)
  {
    cerr << "read_bits_profiler:header record mismatch\n";
    cerr << "Expected '" << HEADER_RECORD << "' got '" << buffer << "'\n";
    return 0;
  }

  number_fixed = 0;
  number_sparse = 0;

  while (input.next_record (buffer))
  {
    if (0 == buffer.length() || buffer.starts_with('#'))
      continue;

    if (buffer.starts_with(COUNT_FIXED))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(number_fixed) || number_fixed < 1)
      {
        cerr << "The number of fixed fingerprints must be a whole +ve number '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(COUNT_SPARSE))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(number_sparse) || number_sparse < 1)
      {
        cerr << "The number of sparse fingerprints must be a whole +ve number '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(BIT_COUNT_CUTOFF))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(bit_count_cutoff) || bit_count_cutoff < 1)
      {
        cerr << "The number of bit count cutoff must be a whole +ve number '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(INCLUDE_LOWER_COUNTS))
    {
      buffer.remove_leading_words(1);
      if (! buffer.numeric_value(include_lower_counts) || include_lower_counts  < 0 )
      {
        cerr << "The include_lower_counts must be 0 or 1 '" << buffer << "'\n";
        return 0;
      }
    }
    else if (buffer.starts_with(CLASS_NAMES))
    {
      buffer.remove_leading_words(1);
      const_IWSubstring token;
      int i=0;
      while( buffer.nextword(token,i))
      {
         class_names.push_back(token);
         class_name_to_id[token]=i;
      }
      number_classes = class_names.size();
    }
    else if ('|' == buffer)   // end of preamble
      break;
    else
    {
      cerr << "Unrecognised record in bits profiler cross reference '" << buffer << "'\n";
      return 0;
    }
  }
  int fpi;
  for(fpi =0;fpi<(number_fixed+number_sparse);++fpi)
  {
    while (input.next_record(buffer))
    {
      if (buffer.starts_with('#') || 0 == buffer.length())
        continue;
        
      if ('|' == buffer)
        break;
        
      if (! read_bits_profiler_record (buffer, bits_profiler,fpi))
      {
        cerr << "Invalid record in bits profiler '" << buffer << "', line " << input.lines_read() << endl;
        return 0;
      }
    }
  }

  return fpi;
}

static int
read_bits_profiler( Bits_Profiler & bits_profiler,const const_IWSubstring & fname)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "GFP_to_Feature_Map::read_existing_mapping:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_bits_profiler (bits_profiler,input);
}

static int calc_sum_frequencies (const Frequencies & freq)
{
   int sum =0;
   for (uint32_t i=0;i<freq.size();i++)
   {
     sum += freq[i];
   }
   return sum;
}
static int
write_bits_profiler (const Bits_Profiler & bits_profiler,
                     IWString_and_File_Descriptor & output)
{
  //int bit_count = bits_profiler.size();
  output << HEADER_RECORD << '\n';
  output << '\n';
  if (number_fixed)
    output << COUNT_FIXED << ' ' << number_fixed << '\n';
  if (number_sparse)
    output << COUNT_SPARSE << ' ' << number_sparse << '\n';
  
  output << CLASS_NAMES;
  for(int i=0;i<number_classes;i++)
    output << ' ' << class_names[i];
  output << '\n';
  
  output << BIT_COUNT_CUTOFF << ' ' << bit_count_cutoff << '\n';
  output << INCLUDE_LOWER_COUNTS << ' ' << include_lower_counts << '\n';
	
  double sum = 0;
  for(uint32_t i=0;i<training_set.size();i++)
    sum += training_set[i];
  for(uint32_t i=0;i<training_set.size();i++)
  {
    training_set[i] /= sum;
  }
  
  int current_fp=-1;

  for (Bits_Profiler::const_iterator it = bits_profiler.begin(); it != bits_profiler.end(); ++it)
  {
    if ((*it).first.first!=current_fp)
    {
      current_fp++;
      output << "|\n";
    }   

    const Frequencies & v = (*it).second;

    double sum = calc_sum_frequencies(v);
    if (min_support_level > 1 && sum < min_support_level)
      continue;

    output << (*it).first.second.first <<' '<<(*it).first.second.second;
    for(int i=0;i<number_classes;i++)
    {
      double a = v[i];
      a = log( ( a + laplace_correction_factor) /( sum * training_set[i] + laplace_correction_factor ) ) ;
      output << ' ' << a; 
    }
    output<< '\n';
    
    output.write_if_buffer_holds_more_than(32768);
  }

  output << "|\n";

  return 1;
}


static int
increment_bit (Bits_Profiler & bits_profiler,
               const int &fpi,
               const unsigned int & b,
               const Class_IDs &act_value,const int c)
{
  FP_And_Bit_And_Count fabc(fpi,Bit_And_Count(b,c));
  
  Bits_Profiler::iterator f = bits_profiler.find(fabc);

  if (f == bits_profiler.end())
  {
    Frequencies freq(number_classes);
    for(uint32_t i=0;i<act_value.size();i++)
      freq[act_value[i]] = 1;
    bits_profiler[fabc] = freq;
  }
  else
  {
    for(uint32_t i=0;i<act_value.size();i++)
      (*f).second[act_value[i]] += 1;
  }
  return 1;
}

static int
profile_bits_set (const IW_General_Fingerprint & fp,
                  Bits_Profiler & bits_profiler,
                  const Class_IDs &act_value)
{
  number_fixed=number_fingerprints();
  number_sparse=number_sparse_fingerprints();
  
  int fpi=0;
  
  while (fpi<number_fixed)
  {
    const IWDYFP & dyfp = fp[fpi];
    int b;
    int i = 0;

    while ((b = dyfp.next_on_bit(i)) >= 0)
    {
      increment_bit(bits_profiler,fpi,b,act_value,1);
    }
    fpi++;
  }
  
  while (fpi<(number_fixed+number_sparse))
  {
    const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(fpi-number_fixed);
          
    int i = 0;
    unsigned int b;
    int c;
  
    while (sfp.next_bit_set(i, b, c))
    {
       if (c > bit_count_cutoff)
         c = bit_count_cutoff;
       while( c >= 1){
         increment_bit(bits_profiler, fpi, b,act_value,c);
         c--;
         if(!include_lower_counts)
           break;
       }
    }
    fpi++;
  }

  return 1;
}            

static int
profile_bits_set (const IW_General_Fingerprint & fp,
                  Bits_Profiler & bits_profiler)
{
  const IWString & id = fp.id();
        
  Class_IDs act_value;
   
  if (! find_identifier_in_activity_hash(activity, id,act_value))
  {
    cerr << "No activity data for '" << id << "'\n";
    if (ignore_missing_activity)
      return 1;

    return 0;
  }
  for(uint32_t i=0;i<act_value.size();i++)
    training_set[act_value[i]] += 1;
  return profile_bits_set(fp,bits_profiler,act_value);
}

static int
profile_bits_set (IW_TDT & tdt,
                  Bits_Profiler & bits_profiler)
{
  IW_General_Fingerprint fp;

  int fatal;

  if (! fp.construct_from_tdt(tdt, fatal))
  {
    cerr << "Cannot build fingerprint\n";
    return 0;
  }

  const_IWSubstring tmp = fp.id();

  if (tmp.nwords() > 1)
  {
    tmp.truncate_at_first(' ');
    fp.set_id(tmp);
  }

//cerr << "Profiling '" << tmp << "'\n";

  return profile_bits_set(fp, bits_profiler);
}

static int
profile_bits_set (iwstring_data_source & input,
                  Bits_Profiler & bits_profiler)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! profile_bits_set (tdt, bits_profiler))
    {
      cerr << "Fatal error processing '" << tdt << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
profile_bits_set (const char * fname,
                  Bits_Profiler & bits_profiler)
{
  iwstring_data_source input (fname);

  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return profile_bits_set (input, bits_profiler);
}

//ref for binary : http://pubs.acs.org/doi/abs/10.1021/ci300435j 


static int
calc_bit_index(const Bits_Profiler & bits_profiler,
               const int &fpi,
               const unsigned int & b, 
               const int c,
               Frequencies &index)
{
  FP_And_Bit_And_Count fabc(fpi,Bit_And_Count(b,c));
  
  Bits_Profiler::const_iterator f = bits_profiler.find(fabc);

  if (f == bits_profiler.end())
    return 0;

  index=(*f).second;
  
  return 1;     
}

static int
bayesian_predict (const IW_General_Fingerprint & fp,
                  const Bits_Profiler & bits_profiler,
                  IWString_and_File_Descriptor &output)
{
  if (number_fixed!=number_fingerprints())
    return 0;
  if (number_sparse!=number_sparse_fingerprints())
    return 0;
        
  const IWString & id = fp.id();
  int fpi=0;
  Frequencies sum(number_classes);
  Frequencies index;
  int n=0;
        
  while (fpi<number_fixed)
  {
    const IWDYFP & dyfp = fp[fpi];
    int b;
    int i = 0;

    while ((b = dyfp.next_on_bit(i)) >= 0)
    {
      if (calc_bit_index(bits_profiler, fpi, b,1, index))
      {
        for(int nc =0;nc<number_classes;nc++)
          sum[nc]+=index[nc];
        ++n;
      }
    }
    fpi++;
  }
        
  while (fpi<(number_fixed+number_sparse))
  {
    const Sparse_Fingerprint & sfp = fp.sparse_fingerprint(fpi-number_fixed);
                
    int i = 0;
    unsigned int b;
    int c;
        
    while (sfp.next_bit_set(i, b, c))
    {
      if (c > bit_count_cutoff)
        c = bit_count_cutoff;
      while( c >= 1){
        if (calc_bit_index(bits_profiler, fpi,b,c,index))
        {
          for(int nc =0;nc<number_classes;nc++)
          {
             sum[nc]+=index[nc];
          }
          ++n;
        }
        c--;
        if(!include_lower_counts)
          break;
      }
    }
    fpi++;
  }

//  sum+=log(training_set_active/training_set_inactive);

//sum+=log((1 + training_set_active)/(1 + training_set_inactive));
  Frequencies prob(number_classes);
  double max_prob = sum[0];
  int max_class = 0;
  for(int i=0;i<number_classes;i++)
  {
    if (sum[i] >max_prob)
    {
      max_prob = sum[i];
      max_class = i;
    }
    if (! write_probabilities_in_01_range)
      continue;
    prob[i] = 0.0;
    for(int j =0;j<number_classes;j++)
      prob[i] += exp(sum[j] -sum[i]); 
    prob[i] = 1.0/prob[i];
    if(prob[i] < 0.001)
      prob[i] = 0; 
  }
  

  output << id << ' ' << class_names[max_class];
  for(int i=0;i<number_classes;i++)
  {
    if(write_probabilities_in_01_range)
      output << ' ' << prob[i] ;
    else
      output << ' ' << sum[i] ;
  }
  output << '\n'; 

  output.write_if_buffer_holds_more_than(32768);
  
  return 1;
}

static int
bayesian_predict (IW_TDT & tdt,
        Bits_Profiler & bits_profiler,
        IWString_and_File_Descriptor &output)
{
  IW_General_Fingerprint fp;

  int fatal;

  if (! fp.construct_from_tdt(tdt, fatal))
  {
    cerr << "Cannot build fingerprint\n";
    return 0;
  }

  const_IWSubstring tmp = fp.id();

  if (tmp.nwords() > 1)
  {
    tmp.truncate_at_first(' ');
    fp.set_id(tmp);
  }

//cerr << "Profiling '" << tmp << "'\n";
  return bayesian_predict (fp, bits_profiler, output);
}

static int
bayesian_predict (iwstring_data_source & input,
                  Bits_Profiler & bits_profiler,
                  IWString_and_File_Descriptor &output)
{
  IW_TDT tdt;
  while (tdt.next(input))
  {
    if (! bayesian_predict (tdt, bits_profiler,output))
    {
      cerr << "Fatal error processing '" << tdt << "'\n";
      return 0;
    }
  }

  return 1;
}


static int
bayesian_predict (const char *fname,
                  Bits_Profiler & bits_profiler,
                  IWString_and_File_Descriptor &output)
{
  iwstring_data_source input(fname);
        
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return bayesian_predict(input, bits_profiler, output);
}       

static void
usage (int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Naive Bayesian training and prediction\n";
  cerr << " -A <fname>     activity file for training\n";
  cerr << " -h             activity file has header line\n";
  cerr << " -B <fname>     bits profiler file for predicting\n";
  cerr << " -c <col>       activity in column <col> of activity file\n";
  cerr << " -F <tag>       tag to process - standard gfp_* syntax\n";
  cerr << " -t <number>    cutoff for bit count (default 1)\n";
  cerr << " -e             do not consider lower bit counts\n";
  cerr << " -l <float>     Laplace correction factor (default 1.0)\n";
  cerr << " -p <number>    min support level for output\n";
  cerr << " -j             write probabilities in 0 to 1 range\n";
  cerr << " -v             verbose output\n";
  
  exit (rc);
}

static int
gfp_naive_bayesian (int argc, char ** argv)
{
  Command_Line cl (argc, argv, "vA:B:c:F:l:p:t:ezjh");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage (1);
  }

  verbose = cl.option_count ('v');

  if (need_to_call_initialise_fingerprints (cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise fingerprint option(s)\n";
      usage(4);
    }
  }
  
  if ( ! cl.option_present('A') && ! cl.option_present('B'))
  {
    cerr << "Must specify activity file via the -A option, or specify bit profiler file via -B option\n";
    usage(4);
  }
  
  if (cl.option_present('c'))
  {
    if (! cl.value('c', activity_column) || activity_column < 1)
    {
      cerr << "The activity column (-c) must be a whole +ve number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Activity data in column " << activity_column << " of activity file\n";

    activity_column--;
  }

  if (cl.option_present('j'))
  {
    write_probabilities_in_01_range = 1;

    if (verbose)
      cerr << "Probabilities written in -1 to 1 range\n";
  }
  if (cl.option_present('h'))
  {
    header_in_activity_file  = 1;

    if (verbose)
      cerr << "Activity file has header line\n";
  }

  if (cl.option_present('z'))
  {
    strip_leading_zeros = 1;

    if (verbose)
      cerr << "Will strip leading zero's from identifiers\n";
  }
  
  if (cl.option_present('A'))
  {
    const_IWSubstring a = cl.string_value('A');

    if (! read_activity_data(a, activity))
    {
      cerr << "Cannot read activity data from '" << a << "'\n";
      return 5;
    }
                
    if (verbose)
      cerr << "Read " << activity.size() << " activity values from '" << a << "'\n";
  }
  
  if (cl.option_present('l'))
  {
    if (! cl.value('l', laplace_correction_factor) || laplace_correction_factor <=0)
    {
      cerr << "The Laplace correction factor (-l) must be a + float number\n";
      usage(4);
    }

    if (verbose)
      cerr << "Laplace correction factor set as " << laplace_correction_factor << " \n";
  }
  
  if (cl.option_present('p'))
  {
    if (! cl.value('p', min_support_level) || min_support_level < 1)
    {
      cerr << "The minimum support level (-p) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will drop bits present in " << min_support_level << " molecules or fewer\n";
  }

  train_mode=1;

  if (cl.option_present('B'))
  {
    const_IWSubstring b = cl.string_value('B');

    if (! read_bits_profiler(bits_profiler,b))
    {
      cerr << "Cannot read bits profiler from '" << b << "'\n";
      return 5;
    }
    
    train_mode=0;
                
    if (verbose)
      cerr << "Read " << bits_profiler.size() << " bits from '" << b << "'\n";
  }

  if (cl.option_present('t'))
  {
    if(0 == train_mode)
    {
      cerr << "Do not support -p option for prediction \n";
      usage(3);
    }
    if (! cl.value('t', bit_count_cutoff) || bit_count_cutoff < 1)
    {
      cerr << "The bit count cutoff (-t) must be a whole +ve number\n";
      usage(3);
    }

    if (verbose)
      cerr << "Will set bit count cutoff to " << bit_count_cutoff << "\n";
  }

  if (cl.option_present('e'))
  {
    if(0 == train_mode)
    {
      cerr << "Do not support -e option for prediction \n";
      usage(3);
    }
    include_lower_counts = 0;

    if (verbose)
      cerr << "Will not include lower bit counts\n";
  }

  set_sparsefp_warn_empty_data(0);

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage (2);
  }

  if (train_mode)
  {
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! profile_bits_set(cl[i], bits_profiler))
      {
        cerr << "Fatal error processing '" << cl[i] << "'\n";
        exit (i + 1);
      }
    }

    if (verbose)
      cerr << "Found " << bits_profiler.size() << " fingerprints\n";

    IWString_and_File_Descriptor output(1);
    write_bits_profiler(bits_profiler,output);
    output.flush();
  }
  else
  {
    IWString_and_File_Descriptor output(1);
    
    output << "ID Prediction";
    for(int i=0;i<number_classes;i++)
      output << ' '<<class_names[i];
    output << '\n';

    set_default_iwstring_double_concatenation_precision(3);
        
    for (int i = 0; i < cl.number_elements(); i++)
    {
      if (! bayesian_predict(cl[i], bits_profiler, output))
      {
        cerr << "Fatal error processing '" << cl[i] << "'\n";
        exit (i + 1);
      }
    }
        
    output.flush();
    
  }
 
  return 0;
}


int
main (int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_naive_bayesian (argc, argv);
  return rc;
}
