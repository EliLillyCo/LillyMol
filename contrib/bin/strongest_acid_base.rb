#!/usr/bin/env ruby

# Example the contents of a file from QupKake that contains pKa values.
# Each molecule might be included multiple times.
# Where there are duplicates reduce to the most acidic and most basic.

class MolData
  def initialize()
    @smi = ""
    @idx = -1
    @acid = 99.0
    @base = -99.0
  end
  def extra(smiles, idx, acid_base, pka)
    if acid_base == 'acidic'
      if pka < @acid
        @smi = smiles
        @acid = pka
        @idx =  idx
      end
    elsif acid_base == 'basic'
      if pka > @base
        @smi = smiles
        @base = pka 
        @idx = idx
      end
    else
      $stderr << "Unrecognised acid/base form #{acid_base}\n"
    end
  end

  def write(id)
    if @acid < 14.0
      $stdout << "#{@smi} #{id} idx:#{@idx} acid_base:acidic pka:#{@acid}\n"
    end
    if @base > 1.0
      $stdout << "#{@smi} #{id} idx:#{@idx} acid_base:basic pka:#{@base}\n"
    end
    
  end
end

def main
  id_to_data = {}
  ARGF.each do |line|
    f = line.chomp.split
    id = f[1]
    idx = f[2].gsub(/idx:/, "").to_i
    acid_base = f[3].gsub(/pka_type:/, "")
    pka = f[4].gsub(/pka:/, "").to_f
    if id_to_data.key?(id)
      id_to_data[id].extra(f[0], idx, acid_base, pka)
    else
      id_to_data[id] = MolData.new()
      id_to_data[id].extra(f[0], idx, acid_base, pka)
    end
  end

  id_to_data.each do |id, data|
    data.write(id)
  end
end

main

