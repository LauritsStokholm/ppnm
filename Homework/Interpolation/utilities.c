/* Compilation of utility functions general for each task */

int binarysearch (int list_size, double* list, double val)
{
  /* ............................................................
  * This binarysearch will given a list of size n
  * x = [x0, x1, x2, ... xn-1]
  * and a value z, return the index value i such that
  * x(i) <= z < x(i+1)
  * ............................................................ */

  // Given an ordered list of n>1 elements, and value (continuous) is in 
  // [x(0), x(n-1)], although not necesarrily contained in list (discrete).
  assert ((list_size > 1) && (val >= list[0]) && (val <= list[list_size-1]));

  // Shrink interval [x(idx_min), x(idx_max)] such value is inside
  int idx_min = 0, idx_max = list_size-1;
  while (idx_max - idx_min > 1)
  {
    // Find midpoint and determine if value is on left or right side
    // [x(i), x(mid)] or [x(mid), x(i+1)]
    int idx_mid = (idx_max + idx_min)/2;
    if ( val > list[idx_mid]) idx_min = idx_mid; else idx_max = idx_mid;
  }
  return idx_min;
}

