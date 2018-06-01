/* ----------------------------------------------------------------------
   
|                                                                   		     
 
   | Swap the bytes of a 4 byte word.                                 
   | Example _byte_swap(0x12345678) = 0x78563412                        
   |                                                                    
   | To do this for short int's:                                        
   |   tmp = (x & 0x00FF);                                             
   |   tmp = ((x & 0xFF00) >> 0x08) | (tmp << 0x08);                    
   ---------------------------------------------------------------------- 
*/

void
rick_higgs_byte_swap (int nw, unsigned int * b)
{
   unsigned int tmp;

   for (int i = 0; i < nw; i++)
   {
     tmp = (b[i] & 0x000000FF);
     tmp = ((b[i] & 0x0000FF00) >> 0x08) | (tmp << 0x08);
     tmp = ((b[i] & 0x00FF0000) >> 0x10) | (tmp << 0x08);
     tmp = ((b[i] & 0xFF000000) >> 0x18) | (tmp << 0x08);

     b[i] = tmp;
   }

   return;
}

void
rick_higgs_byte_swap (unsigned int & b)
{
   unsigned int tmp;

   tmp = (b & 0x000000FF);
   tmp = ((b & 0x0000FF00) >> 0x08) | (tmp << 0x08);
   tmp = ((b & 0x00FF0000) >> 0x10) | (tmp << 0x08);
   tmp = ((b & 0xFF000000) >> 0x18) | (tmp << 0x08);

   b = tmp;

   return;
}
