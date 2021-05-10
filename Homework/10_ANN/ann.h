#ifndef HAVE_ANN_H
#define HAVE_ANN_H

// Define structure for Artificial Neural Network (ANN)
typedef struct
{
  int nneurons; // Number of neurons
  double (*f) (double x, void* params); // Activation function
  vector* parameters; // Neural network parameters (3 per neuron)
} ann;

// For supervised learning
typedef struct
{
  ann* ann;     // The artificial neural network
  vector* x;    // Domain (x-values) for learning
  vector* y;    // Co-domain (y-values) for learning
} my_struct_params;

ann*   ann_alloc    (int  n, double (*f)(double x, void* params) );
void   ann_free     (ann* network);
void   ann_init_p   (ann* network, int xmin, int xmax);
double ann_response (ann* network, double x);

void ann_supervised_train ( ann* network, vector* xdata, vector* ydata,
    double (*cost) (vector* x, void* params), double eps);

#endif
