/* by Adel Javanmard, Andrea Montanari and Federico Ricci-Tersenghi */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#define EPS1 (1e-4)  // convergence tolerance
#define EPS2 (1e-10) // eigensystem tolerance
#define MAX_LINE 10000  // maximum line lenght in graph file

// macros for pseudo-random numbers
#define INV_RAND_MAX (1./RAND_MAX)
#define FRANDOM (INV_RAND_MAX * random())
#define pm1 (FRANDOM < 0.5 ? 1 : -1)

typedef struct {
  unsigned int deg, *neigh;
  double *s;
} nodeStruct;

typedef struct {
  unsigned int numNodes, numEdges, *perm;
  nodeStruct *node;
} graphStruct;

void readGraph(char * filename, graphStruct * pg) {
  int i, n1, n2;
  unsigned int *nodesRead, numReads, size=1<<20;
  char line[MAX_LINE];
  FILE * graphFile = fopen(filename, "r");

  nodesRead = (unsigned int *)calloc(size, sizeof(unsigned int));
  numReads = 0;
  pg->numNodes = 0;
  pg->numEdges = 0;
  while (fgets(line, MAX_LINE, graphFile) != NULL) {
    i = 0;
    while (isspace(line[i])) i++;
    if (line[i] != '#' && line[i] != '\0') {
      if (numReads == size) {
	size *= 2;
	nodesRead = (unsigned int *)realloc(nodesRead, size * sizeof(unsigned int));
      }
      if (sscanf(line + i, "%u %u", nodesRead+numReads, nodesRead+numReads+1) != 2) {
	fprintf(stderr,"\n Error reading graph file on line: %s\n", line);
	exit(EXIT_FAILURE);
      }
      if (nodesRead[numReads] > pg->numNodes) pg->numNodes = nodesRead[numReads];
      numReads++;
      if (nodesRead[numReads] > pg->numNodes) pg->numNodes = nodesRead[numReads];
      numReads++;
      pg->numEdges++;
    }
  }
  pg->numNodes++;
  pg->perm = (unsigned int *)calloc(pg->numNodes, sizeof(unsigned int));
  pg->node = (nodeStruct *)calloc(pg->numNodes, sizeof(nodeStruct));
  for (i = 0; i < pg->numNodes; i++) pg->node[i].deg = 0;
  for (i = 0; i < numReads; i++)
    pg->node[nodesRead[i]].deg++;
  for (i = 0; i < pg->numNodes; i++) {
    pg->node[i].neigh = (unsigned int *)calloc(pg->node[i].deg, sizeof(unsigned int));
    pg->node[i].deg = 0;
  }
  for (i = 0; i < numReads; i += 2) {
    n1 = nodesRead[i];
    n2 = nodesRead[i+1];
    pg->node[n1].neigh[pg->node[n1].deg] = n2;
    pg->node[n2].neigh[pg->node[n2].deg] = n1;
    pg->node[n1].deg++;
    pg->node[n2].deg++;
  }
  fclose(graphFile);
}

double norm(double * v, int m) {
  int i;
  double res = 0.0;
  
  for (i = 0; i < m; i++)
    res += v[i] * v[i];

  return sqrt(res);
}

double scalarProd(double * v1, double * v2, int m) {
  int i;
  double res = 0.0;
  
  for (i = 0; i < m; i++)
    res += v1[i] * v2[i];

  return res;
}

void initSpin(graphStruct * pg, int m) {
  int i, j;
  double tmp;

  for (i = 0; i < pg->numNodes; i++) {
    pg->node[i].s = (double *)calloc(m, sizeof(double));
    if (pg->node[i].deg) {
      for (j = 0; j < m; j++)
	pg->node[i].s[j] = pm1;
      tmp = 1.0 / norm(pg->node[i].s, m);
      for (j = 0; j < m; j++)
	pg->node[i].s[j] *= tmp;
    } else {
      for (j = 0; j < m; j++)
	pg->node[i].s[j] = 0.0;
    }
  }
}

double oneStep(graphStruct * pg, int m, double gamma,
	       double * totalMag, double * localField) {
  int i, j, k, ran, num;
  double tmp, diff, maxDiff, new;
  nodeStruct workNode;
  
  for (i = 0; i < pg->numNodes; i++)
    pg->perm[i] = i;
  num = pg->numNodes;
  for (j = 0; j < m; j++)
    totalMag[j] = 0.0;
  for (i = 0; i < pg->numNodes; i++)
    for (j = 0; j < m; j++)
      totalMag[j] += pg->node[i].s[j];
  maxDiff = 0.0;
  while (num) {
    ran = (int)(FRANDOM * num);
    workNode = pg->node[pg->perm[ran]];
    pg->perm[ran] = pg->perm[--num];
    if (workNode.deg) {
      for (j = 0; j < m; j++)
	localField[j] = -gamma * totalMag[j];
      for (k = 0; k < workNode.deg; k++)
	for (j = 0; j < m; j++)
	  localField[j] += pg->node[workNode.neigh[k]].s[j];
      tmp = norm(localField, m);
      if (tmp != 0.0) {
	tmp = 1.0 / tmp;
	diff = 0.0;
	for (j = 0; j < m; j++) {
	  new = tmp * localField[j];
	  diff += (new - workNode.s[j]) * (new - workNode.s[j]);
	  totalMag[j] += 2.0 * (new - workNode.s[j]);
	  workNode.s[j] = new;
	}
	if (diff > maxDiff) maxDiff = diff;
      }
    }
  }
  return maxDiff;
}

void computeProjections(graphStruct * pg, int m) {
  int i, j, k, num, numProj;
  double **cov, *val, **vec, *new;
  double dist, tmp, sumVal;

  // allocate memory for empirical covariance and its eigensystem
  val = (double *)calloc(m, sizeof(double));
  new = (double *)calloc(m, sizeof(double));
  cov = (double **)calloc(m, sizeof(double *));
  vec = (double **)calloc(m, sizeof(double *));
  for (i = 0; i < m; i++) {
    cov[i] = (double *)calloc(m, sizeof(double));
    vec[i] = (double *)calloc(m, sizeof(double));
  }
  // fill the empirical covariance matrix
  num = 0;
  for (i = 0; i < pg->numNodes; i++)
    if (pg->node[i].deg) {
      num++;
      for (j = 0; j < m; j++)
	for (k = 0; k < m; k++)
	  cov[j][k] += pg->node[i].s[j] * pg->node[i].s[k];
    }
  for (j = 0; j < m; j++)
    for (k = 0; k < m; k++)
      cov[j][k] /= num;
  // compute the eiegnsystem by matrix powers
  i = 0;
  sumVal = 0.0;
  while (i < m && sumVal < 0.999999) {
    for (j = 0; j < m; j++)
      vec[i][j] = pm1 / sqrt(m);
    do {
      for (j = 0; j < m; j++)
	new[j] = scalarProd(cov[j], vec[i], m);
      // ortogonalize w.r.t. previous eigenvectors
      for (k = 0; k < i; k++) {
	tmp = scalarProd(new, vec[k], m);
	for (j = 0; j < m; j++)
	  new[j] -= tmp * vec[k][j];
      }
      val[i] = scalarProd(new, vec[i], m);
      tmp = 1.0 / norm(new, m);
      dist = 0.0;
      for (j = 0; j < m; j++) {
	new[j] *= tmp;
	dist += (new[j] - vec[i][j]) * (new[j] - vec[i][j]);
	vec[i][j] = new[j];
      }
    } while (dist > EPS2);
    sumVal += val[i];
    i++;
  }
  numProj = i;
  printf("\n\n# Eigenvalues of the empirical covariance:");
  for (i = 0; i < numProj; i++)
    printf(" %g", val[i]);
  printf("\n#\n");
  printf("# nodeIndex, degree, projections along the first %i eigenvectors\n", numProj);
  for (i = 0; i < pg->numNodes; i++) {
    if (pg->node[i].deg) {
      printf("%i %i", i, pg->node[i].deg);
      for (j = 0; j < numProj; j++)
	printf(" %g", scalarProd(pg->node[i].s, vec[j], m));
      printf("\n");
    } else {
      printf("%i 0\n", i);
    }
  }
}

int main(int argc, char *argv[]) {
  unsigned int m, maxIter, randomSeed, trueNumNodes, i, t;
  double gamma, maxDiff, *totalMag, *localField;
  graphStruct graph;
  
  if (argc != 5 && argc != 6) {
    fprintf(stderr, "usage: %s <graphFile> <m> <maxIter> <randomSeed> [<gamma>]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  printf("# Reading graph file %s...", argv[1]); fflush(stdout);
  readGraph(argv[1], &graph);
  printf("done.\n");
  m = atoi(argv[2]);
  maxIter = atoi(argv[3]);
  randomSeed = atoi(argv[4]);
  trueNumNodes = 0;
  for (i = 0; i < graph.numNodes; i++) if (graph.node[i].deg) trueNumNodes++;
  if (argc == 6) gamma = atof(argv[5]);
  else gamma = 2.0 * graph.numEdges / trueNumNodes / trueNumNodes;

  printf("# running SDP clustering on data in file %s\n", argv[1]);
  printf("# num nodes = %i   num edges = %i\n", trueNumNodes, graph.numEdges);
  printf("# m = %i   maxIter = %i   gamma = %g   randomSeed = %i\n",
	 m, maxIter, gamma, randomSeed);
  printf("# EPS1 = %g (convergence tolerance)  EPS2 = %g (eigensystem tolerance)\n", EPS1, EPS2);
  srandom(randomSeed);
  totalMag = (double *)calloc(m, sizeof(double));
  localField = (double *)calloc(m, sizeof(double));
  initSpin(&graph, m);
  printf("#\n# time, maxDiff\n");
  t = 0;
  do {
    maxDiff = oneStep(&graph, m, gamma, totalMag, localField);
    t++;
    printf("%i %g\n", t, maxDiff);
  } while (t < maxIter && maxDiff > EPS1);
  computeProjections(&graph, m);
  return EXIT_SUCCESS;
}
