#include <math.h>
#include <stdio.h>

#define PARTICION 11  // Tamaño de la partición
#define N 9           // Número de incógnitas
#define VISUALIZA 1   // (0) No visualiza la salida, otro valor la visualiza


/*
 * Este programa resuelve numéricamente la ecuación diferencial parcial
 * que describe un problema físico. En este caso, estamos resolviendo una 
 * ecuación de frontera de segundo orden en una dimensión, usando métodos
 * de diferencia finita (específicamente el método de Jacobi y Gauss-Seidel)
 * para aproximar la solución de la ecuación de Poisson:
 *
 *   d²U(x) / dx² = f(x)
 *
 * con condiciones de frontera en los extremos del dominio [xi, xf]:
 * 
 *   U(xi) = vi  y   U(xf) = vf.
 * 
 * La solución numérica se obtiene al discretizar el dominio en una malla de
 * tamaño `h` y resolver el sistema de ecuaciones lineales resultante.
 *
 * Los resultados son guardados en un archivo "resultados.dat" con los valores
 * de la solución `U(x)` para cada valor de `x` en el dominio discretizado.
 * 
 * El archivo de salida tiene dos columnas: la primera es el valor de `x`
 * y la segunda es el valor de la solución `U(x)` calculado en ese punto.
 */


// Lado derecho de la ecuación diferencial parcial
double LadoDerecho(double x)
{
   double pi = 3.1415926535897932384626433832;
   return -pi * pi * cos(pi * x);  // Definición de la función LadoDerecho
}

// Resuelve Ax=b usando el método Jacobi
void Jacobi(double A[N][N], double x[], double b[], int n, int iter)
{
   int i, j, m;
   double sum;
   double xt[N];  // Vector temporal para almacenar los resultados intermedios

   for (m = 0; m < iter; m++)  // Iteraciones del método Jacobi
   {
      for (i = 0; i < n; i++)  // Recorrer cada incógnita
      {
         sum = 0.0;
         for (j = 0; j < n; j++)  // Sumar los términos de la ecuación
         {
            if (i == j) continue;  // Ignorar el término de la diagonal principal
            sum += A[i][j] * x[j];  // Sumar el producto de la matriz y el vector
         }
         if (A[i][i] == 0.0) return;  // Si la diagonal principal es cero, terminamos
         xt[i] = (1.0 / A[i][i]) * (b[i] - sum);  // Resolver para el valor de x[i]
      }
      for (i = 0; i < n; i++)  // Actualizar el vector x
         x[i] = xt[i];
   }
}

// Resuelve Ax=b usando el método Gauss-Seidel
void Gauss_Siedel(double A[N][N], double x[], double b[], int n, int iter)
{
   int i, j, m;
   double sum;

   for (m = 0; m < iter; m++)  // Iteraciones del método Gauss-Seidel
   {
      for (i = 0; i < n; i++)  // Recorrer cada incógnita
      {
         sum = 0.0;
         for (j = 0; j < n; j++)  // Sumar los términos de la ecuación
         {
            if (i == j) continue;  // Ignorar el término de la diagonal principal
            sum += A[i][j] * x[j];  // Sumar el producto de la matriz y el vector
         }
         if (A[i][i] == 0.0) return;  // Si la diagonal principal es cero, terminamos
         x[i] = (1.0 / A[i][i]) * (b[i] - sum);  // Resolver para el valor de x[i]
      }
   }
}

int main()
{
   double xi = -1.0;               // Inicio del dominio
   double xf = 2.0;                // Fin del dominio
   double vi = -1.0;               // Valor en la frontera xi
   double vf = 1.0;                // Valor en la frontera xf
   int n = PARTICION;              // Partición
   double h = (xf - xi) / (n - 1); // Incremento en la malla
   int i;

   double A[N][N];  // Matriz A
   double b[N];     // Vector b
   double x[N];     // Vector x

   double R = 1 / (h * h);   // Factor de la ecuación diferencial
   double P = -2 / (h * h);  // Factor de la ecuación diferencial
   double Q = 1 / (h * h);   // Factor de la ecuación diferencial

   // Primer renglón de la matriz A y vector b
   A[0][0] = P;
   A[0][1] = Q;
   b[0] = LadoDerecho(xi) - vi * R;
   
   // Renglones intermedios de la matriz A y vector b
   for (i = 1; i < N - 1; i++)
   {
      A[i][i - 1] = R;
      A[i][i] = P;
      A[i][i + 1] = Q;
      b[i] = LadoDerecho(xi + h * i);
   }
   
   // Renglón final de la matriz A y vector b
   A[N - 1][N - 2] = R;
   A[N - 1][N - 1] = P;
   b[N - 1] = LadoDerecho(xi + h * N - 2) - vf * Q;

   // Resuelve el sistema lineal Ax=b utilizando Gauss-Seidel
   Gauss_Siedel(A, x, b, N, 1000);
   
   // Resuelve el sistema lineal Ax=b utilizando Jacobi
   Jacobi(A, x, b, N, 1000);

   // Crear y abrir el archivo para escribir los resultados
   FILE *file = fopen("resultados.dat", "w");  // Abrir archivo en modo escritura
   if (file == NULL)
   {
       printf("Error al abrir el archivo para escritura.\n");
       return 1;  // Terminar el programa si hubo error al abrir el archivo
   }

   // Escribir los resultados de x y U(x) en el archivo
   for (i = 0; i < N; i++)
   {
       fprintf(file, "%f %f\n", xi + i * h, x[i]);  // Escribir en el archivo
   }

   // Cerrar el archivo después de escribir los resultados
   fclose(file);

   printf("Resultados escritos en 'resultados.dat'\n");  // Confirmación de que se guardaron los resultados

   return 0;
}
