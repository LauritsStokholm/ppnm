struct data {vector* energy; vector* cross_section; vector* error;};
typedef struct data data;



double
Breit_Wigner (double E, double m, double G, double A)
{
  //double E = vector_get (x, 0); // Energy
  //double m = vector_get (x, 1); // Mass
  //double G = vector_get (x, 2); // Width of sought resonance
  //double A = vector_get (x, 3); // Scale factor
  return A / ( (E-m)*(E-m) + (G*G/4.));
}


// Fitting the data by minimising the deviation function
double
cost (vector* x, void* params)
{
  double m = vector_get (x, 0);
  double G = vector_get (x, 1);
  double A = vector_get (x, 2);

  data* parameters = (data*) params;
  vector* E      = parameters -> energy;
  vector* sigma  = parameters -> cross_section;
  vector* dsigma = parameters -> error;

  double sum = 0;
  for (int i=0; i<E->size; i++)
  {
    double Ei = vector_get (E, i);
    double si = vector_get (sigma, i);
    double di = vector_get (dsigma, i);

    sum += pow(Breit_Wigner (Ei, m, G, A) - si, 2) / (di*di);
  }
  return sum;
}


void
task_b (void)
{

  // Read datafile into vector
  FILE* datafile = fopen ("higgs.dat", "r");
  FILE* fp       = fopen ("output.txt", "ab"); // dont overwrite

  // scan data into vector
  double Ei, Si, Di;
  // Prepare size of vector
  int nsteps, dim = 0;

  // Read header (to set pointer beyond)
  fscanf(datafile, "%*[^\n]\n");
  while (fscanf (datafile, "%lg %lg %lg\n", &Ei, &Si, &Di) != EOF) { dim ++; }

  // Load data into vectors
  vector* E = vector_alloc (dim);
  vector* S = vector_alloc (dim);
  vector* D = vector_alloc (dim);

  int i=0;
  rewind (datafile);
  fscanf(datafile, "%*[^\n]\n");
  while (fscanf (datafile, "%lg\t%lg\t%lg\n", &Ei, &Si, &Di) != EOF)
  {
    vector_set (E, i, Ei);
    vector_set (S, i, Si);
    vector_set (D, i, Di);
    i++;
  }

  // Prepare parameters
  data p = {E, S, D};

  // Initial guess (starting point) for vector x = (m, G, A)
  vector* x = vector_alloc (3);

  // Inspired by lecturer
  vector_set (x, 0, 120);
  vector_set (x, 1, 4);
  vector_set (x, 2, 5);

  // epsilon: |grad|<eps on return
  double eps = 1e-12;
  nsteps = qnewton (cost, x, (void*) &p, eps);

  fprintf (fp, "\n\n\nHIGGS BOSON\n");
  fprintf (fp, "Using qnewton\n");
  nsteps = qnewton (cost, x, (void*) &p, eps);
  fprintf_results  (fp, nsteps, x, (void*) &p, cost);
  fprintf (fp, "where the parameters are on the form [mass, witdh of sought resonance, width of the Higgs boson]\n");

  fprintf (fp, "\nUsing nmsimplex\n");
  vector* start = vector_alloc (3);
  vector* step  = vector_alloc (3);

  vector_set (start, 0, 120); vector_set (start, 1, 4); vector_set (start, 2, 5);
  vector_set (step, 0, 1) ; vector_set (step, 1, 0.1); vector_set (step, 2, 0.1);

  double simplex_step_goal = 1e-5;
  nsteps = nmsimplex (cost, start, step, (void*) &p, simplex_step_goal);
  fprintf_results    (fp, nsteps, start, (void*) &p, cost);
  fprintf (fp, "where the parameters are on the form [mass, witdh of sought resonance, width of the Higgs boson]\n");

  vector_free (start);
  vector_free (step);

  vector_free (x);
  vector_free (E);
  vector_free (S);
  vector_free (D);
  fclose (datafile);
  fclose (fp);
  return;
}
