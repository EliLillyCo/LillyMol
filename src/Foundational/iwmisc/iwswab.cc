void
rick_higgs_byte_swap (int nw, unsigned int * bi)
{
  nw = nw * sizeof (int);

  unsigned char * b = reinterpret_cast<unsigned char *> (bi);

  register unsigned char t;

  for (int i = 0; i < nw; i += 4)
  {
    t = *b;
    *b = b[3];
    b[3] = t;
    t = b[1];
    b[1] = b[2];
    b[2] = t;

    b += 4;
  }

  return;
}
