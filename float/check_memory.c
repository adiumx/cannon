#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include "memory_parallel.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
// Modifica esto para el tamaño de matriz
#define N 2048
#define TRUE 1
#define FALSE 0

void imprimematriz(double *matriz, int n);
void MatrixMultiply(int n, double *a, double *b, double *c);
double *create_array_as_matrix(int r, int c);
void populate_array_as_matrix(double *arr, int r, int c);
int array_as_matrix_equals(double *a, double *b, int r, int c);

double *a;
double *b;
double *c;

double bloque_a[16][N * N];
double bloque_b[16][N * N];
double bloque_c[16][N * N];
void MatrixMatrixMultiply(int n, double *a, double *b, double *c, double *c_grande, MPI_Comm comm)
{
    int i;
    int nlocal;
    int npes, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank, coords[2];
    int shiftsource, shiftdest;
    MPI_Status status;
    MPI_Comm comm_2d;

    /* Obtiene toda la información relacionada al canal de comunicación */
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &myrank);

    /* Crea la topología cartesiana */
    dims[0] = dims[1] = sqrt(npes);

    /* Crea los periodos para las conecciones*/
    periods[0] = periods[1] = 1;

    /* Crea la topología cartesiana, reordenando las tareas */
    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

    /* Obtiene las tareas y coordenadas segun la topología*/
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

    /* Aqui se mueven intercambia el contenido de las tareas a la izquierda o arriba
        segun sea la matriz A o B.

        Esto es parte fundamental del algoritmo
     */
    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    /* Determina la dimension de la matriz de bloque local */
    nlocal = n / dims[0];

    /* Inicializa el orden inicial de las matrices A y B*/
    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);

    /* Entra en el loop */
    for (i = 0; i < dims[0]; i++)
    {
        /*Aqui multiplica como lo hace cannon*/
        MatrixMultiply(nlocal, a, b, c); /*c=c+a*b*/
        /* Mueve las columnas de A a la izquierda una casilla */
        MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_2d, &status);
        /* Mueve las filas de B arriba una casilla */
        MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_DOUBLE, uprank, 1, downrank, 1, comm_2d, &status);
    }

    /* Regresa la distribución original de A y B */
    MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_2d, &status);
    /* Limpia el comunicador */
    MPI_Comm_free(&comm_2d);

    int ind;
    int n_grande = n;
    int ind_col;
    for (ind = 0; ind < nlocal; ind++)
    {
        for (ind_col = 0; ind_col < nlocal; ind_col++)
        {
            int fila_grande = mycoords[0] * nlocal + ind;
            int columna_grande = mycoords[1] * nlocal + ind_col;
            c_grande[fila_grande * n_grande + columna_grande] = c[ind * nlocal + ind_col];
        }
    }

    /* Esta parte es la que va relacionando las submatrices con la matriz más grande */

    if (myrank != 0)
    {
        MPI_Reduce(c_grande, c_grande, n_grande * n_grande, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
    else
    {
        MPI_Reduce(MPI_IN_PLACE, c_grande, n_grande * n_grande, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
}

void MatrixMultiply(int n, double *a, double *b, double *c)
{
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

void imprimematriz(double *matriz, int n)
{
    int i;
    int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%0.2f\t", matriz[i * n + j]);
        }
        printf("\n");
    }
}


int main (int argc, char *argv[])
{
	 struct timeval start2, stop;
    a = create_array_as_matrix(N, N);
    b = create_array_as_matrix(N, N);
    c = create_array_as_matrix(N, N);

    populate_array_as_matrix(&a[0], N, N);
    populate_array_as_matrix(&b[0], N, N);
    int my_rank, nro_procesos;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int processor_name_len;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Get_processor_name(processor_name, &processor_name_len);
    int n_chica = N / sqrt(nro_procesos);

    int max_fila_bloque = n_chica;
    int max_columna_bloque = n_chica;

    int fila_bloque;
    int columna_bloque;
    int fila;
    int columna;

    int n_bloques = N / n_chica;
    columna_bloque = my_rank % n_bloques;
    fila_bloque = (my_rank - columna_bloque) / n_bloques;

    int indice_bloque = 0;

    for (fila = fila_bloque * n_chica; fila < fila_bloque * n_chica + n_chica; fila++)
    {
        for (columna = columna_bloque * n_chica; columna < columna_bloque * n_chica + n_chica; columna++)
        {
            bloque_a[my_rank][indice_bloque] = a[fila * N + columna];
            bloque_b[my_rank][indice_bloque] = b[fila * N + columna];
            bloque_c[my_rank][indice_bloque] = 0;
            indice_bloque++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /* Comienza a correr el tiempo */
    gettimeofday(&start2, 0);
    double start = MPI_Wtime();
    /* Aqui se encuentra el algoritmo de cannon com tal */
    MatrixMatrixMultiply(N, bloque_a[my_rank], bloque_b[my_rank], bloque_c[my_rank], &c[0], MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        // Se detiene el tiempo y se muestran los resultados

        /* -----CANNON----- */

        double end = MPI_Wtime();
        // printf("\nTiempo Cannon: %.4f segundos\n", (end - start));
        gettimeofday(&stop, 0);
	double t=0;
	//fprintf(stdout, "%.8f\n",
          //      (stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6));
	t=(stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6);
	FILE* fichero;     
	if(N==16){
		fichero = fopen("ftiempo16.txt", "at");
	}
	if(N==32){
		fichero = fopen("ftiempo32.txt", "at");
		
	}
	if(N==64){
		fichero = fopen("ftiempo64.txt", "at");
		
	}
	if(N==128){
		fichero = fopen("ftiempo128.txt", "at");
		
	}
	if(N==256){
		fichero = fopen("ftiempo256.txt", "at");
	}
		
	if(N==512){
		fichero = fopen("ftiempo512.txt", "at");
		
	}
	if(N==1024){
		fichero = fopen("ftiempo1024.txt", "at");
		
	}
	if(N==2048){
		fichero = fopen("ftiempo2048.txt", "at");
	}
    	fprintf (fichero, "%.8f\n", t);
    	fclose(fichero);

        // /* -----SECUENCIAL------ */
        // /* SOLO PARA COMPROBAR */

        // double *d = create_array_as_matrix(N, N);
        // /* Mide el tiempo del algoritmo secuencial */
        // gettimeofday(&start2, 0);
        // start = MPI_Wtime();
        // MatrixMultiply(N, a, b, d);
        // /* Si queremos podemos imprimirla */
        //    imprimematriz(&d[0], N);
        // printf("Cannon \n\n");
        //    imprimematriz(&c[0], N);

        // int equal = array_as_matrix_equals(&d[0], &c[0], N, N);
        // /* Este punto es importantisimo

        //     Se evalua si efectivamente las dos matrices (Cannon y Secuencial)
        //     estan dando el mismo resultado.

        //  */
        // if (equal)
        // {
        //     /* Si son iguales, quiere decir que cannon esta bien ejecutado

        //         y Muestra el tiempo del Algoritmo secuencial

        //      */
        //     printf("\nSon iguales\n");
        //     end = MPI_Wtime();
        //     printf("\nTiempo Secuencial: %.4f segundos\n", (end - start));
        //     gettimeofday(&stop, 0);
        //     fprintf(stdout, "Time = %.6f\n\n",
        //             (stop.tv_sec + stop.tv_usec * 1e-6) - (start2.tv_sec + start2.tv_usec * 1e-6));
        // }
        // else
        // {
        //     printf("\n No son iguales\n");
        // }
    }
    //printf("Number_of_processes=%03d, My_rank=%03d, processor_name=%5s\n", 
       // nro_procesos, my_rank, processor_name);

    

    long vmrss_per_process[nro_procesos];
    long vmsize_per_process[nro_procesos];
    get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, nro_procesos);
	FILE* fichero1; 
	FILE* fichero2; 
	FILE* fichero3; 
	FILE* fichero4; 
	if(N==16){
		fichero1 = fopen("VmRSS16.txt", "at");
		fichero2 = fopen("VmSize16.txt", "at");
	}
	if(N==32){
		fichero1 = fopen("VmRSS32.txt", "at");
		fichero2 = fopen("VmSize32.txt", "at");
	}
	if(N==64){
		fichero1 = fopen("VmRSS64.txt", "at");
		fichero2 = fopen("VmSize64.txt", "at");
	}
	if(N==128){
		fichero1 = fopen("VmRSS128.txt", "at");
		fichero2 = fopen("VmSize128.txt", "at");
	}
	if(N==256){
		fichero1 = fopen("VmRSS256.txt", "at");
		fichero2 = fopen("VmSize256.txt", "wt");
	}
	if(N==512){
		fichero1 = fopen("VmRSS512.txt", "at");
		fichero2 = fopen("VmSize512.txt", "at");
	}
	if(N==1024){
		fichero1 = fopen("VmRSS1024.txt", "at");
		fichero2 = fopen("VmSize1024.txt", "at");
	}
	if(N==2048){
		fichero1 = fopen("VmRSS2048.txt", "at");
		fichero2 = fopen("VmSize2048.txt", "at");
	}
	
    if (my_rank == 0)
    {
        for (int k = 0; k < nro_procesos; k++)
        {
            //printf("Process %03d: \nVmRSS = %6ld KB\nVmSize = %6ld KB\n", 
              //  k, vmrss_per_process[k], vmsize_per_process[k]);
    		fprintf (fichero1, "%ld\n",vmrss_per_process[k] );
		fprintf (fichero2, "%ld\n",vmsize_per_process[k]);
        }	
    }
	if(N==16){
		fichero3 = fopen("globalvmRSS16.txt", "at");
		fichero4 = fopen("globalvmSize16.txt", "at");
	}
	if(N==32){
		fichero3 = fopen("globalvmRSS32.txt", "at");
		fichero4 = fopen("globalvmSize32.txt", "at");
	}
	if(N==64){
		fichero3 = fopen("globalvmRSS64.txt", "at");
		fichero4 = fopen("globalvmSize64.txt", "at");
	}
	if(N==128){
		fichero3 = fopen("globalvmRSS128.txt", "at");
		fichero4 = fopen("globalvmSize128.txt", "at");
	}
	if(N==256){
		fichero3 = fopen("globalvmRSS256.txt", "at");
		fichero4 = fopen("globalvmSize256.txt", "at");
	}
	if(N==512){
		fichero3 = fopen("globalvmRSS512.txt", "at");
		fichero4 = fopen("globalvmSize512.txt", "at");
	}
	if(N==1024){
		fichero3 = fopen("globalvmRSS1024.txt", "at");
		fichero4 = fopen("globalvmSize1024.txt", "at");
	}
	if(N==2048){
		fichero3 = fopen("globalvmRSS2048.txt", "at");
		fichero4 = fopen("globalvmSize2048.txt", "at");
	}
    long global_vmrss, global_vmsize;
    get_global_memory_usage_kb(&global_vmrss, &global_vmsize, nro_procesos);
    if (my_rank == 0)
    {
        //printf("Global memory usage: VmRSS = %6ld KB\nVmSize = %6ld KB\n", 
          //  global_vmrss, global_vmsize);
		fprintf (fichero3, "%ld\n",global_vmrss );
		fprintf (fichero4, "%ld\n",global_vmsize );
    }
	fclose(fichero1);  
	fclose(fichero2); 
	fclose(fichero3);  
	fclose(fichero4);     
    MPI_Finalize();
    return 0;
}
double *create_array_as_matrix(int r, int c)
{
    double *mat = calloc(r * c, sizeof(double));
    return mat;
}

void populate_array_as_matrix(double *arr, int r, int c)
{
    int j;
    for (j = 0; j < r * c; j++)
    {
        arr[j] = rand() % 2000 + 1000;
    }
}

int array_as_matrix_equals(double *a, double *b, int r, int c)
{
    int i = 0;
    for (i = 0; i < r * c; i++)
    {
        if (a[i] != b[i])
        {
            return FALSE;
        }
    }
    return TRUE;
}
