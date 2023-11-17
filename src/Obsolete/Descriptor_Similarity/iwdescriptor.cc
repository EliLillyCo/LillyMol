
/*
  Sometimes we have numbers in the form 123E+4
*/

static int big_e = 0;

void
set_descriptors_may_contain_big_E (int d)
{
  big_e = d;

  return;
}

int
descriptors_may_contain_big_E() {
  return big_e;
}
