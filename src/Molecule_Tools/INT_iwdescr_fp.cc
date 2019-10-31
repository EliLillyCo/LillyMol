void
Descriptor::produce_fingerprint (int bitnum, Sparse_Fingerprint_Creator & sfc) const
{
  if (0 == _fingerprint_replicates)
  {
    cerr << "Descriptor::produce_fingerprint:not initialised '" << _name << "'\n";
    return;
  }

  float v;

  if (! Set_or_Unset<float>::value(v))
    return;

  int c;

  if (v <= _min)
    c = 1;
  else
  {
    if (v >= _max)
      v = _max;

    c = static_cast<int>((v - _min) / _dy + 0.5f);

    c++;   // ensure non zero
  }
  
//cerr << "Descriptor::produce_fingerprint: descriptor " << _name << " value " << v << " bit " << bitnum << " value " << c << " rep " << _fingerprint_replicates << endl;

  if (1 == _fingerprint_replicates)    // presumably a common case
  {
    sfc.hit_bit(bitnum, c);
    return;
  }

  for (int i = 0; i < _fingerprint_replicates; ++i)
  {
    sfc.hit_bit(bitnum + i, c);
  }

  return;
}
