#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>

#include "iwconfig.h"

#ifdef UNIX
#include <unistd.h>
#endif

#include <assert.h>

using namespace std;

#include "iwrandom.h"

static MTRand_int32 default_random_number_state;

random_number_t
iwrandom ()
{
  return default_random_number_state.closed_open();
}

int
intbtwij (int low, int high)
{
  return low + static_cast<int>(default_random_number_state.closed_closed() * static_cast<double>(high - low + 1));
}

void
iw_set_rnum_seed (random_number_seed_t seed)
{
  default_random_number_state.seed(seed);

  return;
}

static void
jumble_bytes (char * b, int nb,
              int bfrom, int bto)
{
  char bsave = b[bto];

  b[bto] = b[bfrom];
  b[bfrom] = bsave;

  return;
}

random_number_seed_t
random_seed_based_on_time_and_pid()
{
  unsigned int tt = static_cast<unsigned int>(time(0));

  jumble_bytes(reinterpret_cast<char *>(&tt), sizeof(tt), 0, 3);

  int pid = static_cast<int>(IW_GETPID ());

  jumble_bytes(reinterpret_cast<char *>(&pid), sizeof(pid), 2, 3);

  random_number_seed_t seed = static_cast<random_number_seed_t>(tt | pid);

  return seed;
}

random_number_seed_t
random_seed_from_dev_random ()
{
  int randomData = open("/dev/random", O_RDONLY);
  if (randomData < 0)
    return random_seed_based_on_time_and_pid();

  random_number_seed_t rc;
  read(randomData, &rc, sizeof(rc));
  close(randomData);

  return rc;
}

random_number_seed_t
iw_random_seed (void)
{
  random_number_seed_t seed = random_seed_based_on_time_and_pid();

//cerr << "Using seed " << seed << '\n';

  default_random_number_state.seed(seed);

  return seed;
}

Random_Number_Working_Storage::Random_Number_Working_Storage()
{
  MTRand_int32::seed(iw_random_seed());
}

random_number_seed_t
Random_Number_Working_Storage::choose_random_seed()
{
  random_number_seed_t seed = random_seed_based_on_time_and_pid();

  MTRand_int32::seed(seed);

  return seed;
}
