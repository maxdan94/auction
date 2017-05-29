/*
Maximilien Danisch
May 2017
http://bit.ly/maxdan94
maximilien.danisch@gmail.com

Info:
Feel free to use these lines as you wish. This is an efficient C implementation of the Auction algorithm for the assignement problem such as detailed in:
D.P. Bertsekas, A distributed algorithm for the assignment problem,Laboratory for Information and Decision Systems Working Paper (M.I.T., March 1979).

It should easily scale to hundreds of millions of nodes and/or edges if the data is not too adverserial.

To compile:
"gcc auction.c -O3 -o auction".

To execute:
"./auction edgelist.txt res.txt eps".
"edgelist.txt" should contain the edges with integral weights in the bipartite graph: one edge on each line separated by spaces "n1 n2 w", n1 for 0 to n1max and n2 for 0 to n2max.
"res.txt" contains the results: one edge "n1 n2 w" of the assignement on each line.
eps is the stepsize, if eps<1/min(n1,n2) then the algorithm is optimal.
Will print some information in the terminal.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define NLINKS 1000000000

typedef struct {
	unsigned u;//first node
	unsigned v;//second node
	unsigned w;//weight of the edge
} edge;

typedef struct {

	//edge list
	unsigned n1;//number of nodes
	unsigned n2;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	bool b;//b=1 iff n1>n2
	//graph
	unsigned *cd;//cumulative degree: (start with 0) length=dim+1
	unsigned *adj;//list of neighbors
	unsigned *w;//weight of the edges

} bipgraph;

void freebg(bipgraph *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g->w);
	free(g);
}

//compute the maximum of two unsigned
inline unsigned max2u(unsigned a,unsigned b){
	return (a>b) ? a : b;
}

//reading the edgelist of bipartit weighted graph from file
bipgraph* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	bipgraph *g=malloc(sizeof(bipgraph));
	FILE *file;

	g->n1=0;
	g->n2=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u %u\n", &(g->edges[g->e].u), &(g->edges[g->e].v),&(g->edges[g->e].w))==3) {
		g->n1=max2u(g->n1,g->edges[g->e].u);
		g->n2=max2u(g->n2,g->edges[g->e].v);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n1++;
	g->n2++;
	g->edges=realloc(g->edges,g->e*sizeof(edge));
	return g;
}

//Building a special sparse graph structure
void mkbg(bipgraph *g){
	unsigned i,u;
	unsigned *d;

	if (g->n1<=g->n2){
		g->b=0;
		d=calloc(g->n1,sizeof(unsigned));
		g->cd=malloc((g->n1+1)*sizeof(unsigned));
		g->adj=malloc((g->e)*sizeof(unsigned));
		g->w=malloc((g->e)*sizeof(unsigned));
		for (i=0;i<g->e;i++) {
			d[g->edges[i].u]++;
		}
		g->cd[0]=0;
		for (i=1;i<g->n1+1;i++) {
			g->cd[i]=g->cd[i-1]+d[i-1];
			d[i-1]=0;
		}
		for (i=0;i<g->e;i++) {
			u=g->edges[i].u;
			g->adj[g->cd[u] + d[u] ]=g->edges[i].v;
			g->w[g->cd[u] + d[u]++ ]=g->edges[i].w;
		}
	}
	else{
		g->b=1;
		d=calloc(g->n2,sizeof(unsigned));
		g->cd=malloc((g->n2+1)*sizeof(unsigned));
		g->adj=malloc((g->e)*sizeof(unsigned));
		g->w=malloc((g->e)*sizeof(unsigned));
		for (i=0;i<g->e;i++) {
			d[g->edges[i].v]++;
		}
		g->cd[0]=0;
		for (i=1;i<g->n2+1;i++) {
			g->cd[i]=g->cd[i-1]+d[i-1];
			d[i-1]=0;
		}
		for (i=0;i<g->e;i++) {
			u=g->edges[i].v;
			g->adj[g->cd[u] + d[u] ]=g->edges[i].u;
			g->w[g->cd[u] + d[u]++ ]=g->edges[i].w;
		}
		i=g->n1;
		g->n1=g->n2;
		g->n2=i;
	}
	free(d);
	free(g->edges);
}

//bidding towards convergence
unsigned bidding(bipgraph *g,double eps) {
	unsigned u,v,w,i,j,k,res;
	double *p=calloc(g->n2,sizeof(double));//p[i]=price of object i
	unsigned *a=malloc(g->n2*sizeof(unsigned));//a[target]=-1 if not assigned = source if assigned
	unsigned *aw=malloc(g->n2*sizeof(unsigned));//aw[target]=? if not assigned = weigth of edge source-target if assigned
	for (v=0;v<g->n2;v++){
		a[v]=-1;
	}
	//set data structure
	unsigned n_set=g->n1;//
	unsigned *l_set=malloc(g->n1*sizeof(unsigned));//
	//bool *t_set=malloc(g->n1*sizeof(bool));//
	for (u=0;u<n_set;u++){
		l_set[u]=u;
		//t_set[u]=1;
	}

	//double eps=0.1;//1./((double)(g->n1)+1.);///////////
	double max,max2;
	unsigned wmax;

	res=0;
	while (n_set>0){
		u=l_set[n_set-1];
		//printf("n_set, u, w = %u, %u, %u\n",n_set,u,res);
		max=0;
		for (i=g->cd[u];i<g->cd[u+1];i++){
			v=g->adj[i];
			w=g->w[i];
			//printf("w=%u\n",w);
			if ((double)w-p[v]>max){
				max=(double)w-p[v];
				wmax=w;
				k=v;
			}
		}
		//printf("max=%f\n",max);
		if (max==0){
			n_set--;
			continue;
		}
		max2=0;
		for (i=g->cd[u];i<g->cd[u+1];i++){
			v=g->adj[i];
			if (v!=k){
				w=g->w[i];
				if ((double)w-p[v]>max2){
					max2=(double)w-p[v];
				}
			}
		}
		//printf("max2=%f\n",max2);
		p[k]+=max-max2+eps;
		//printf("k=%u\n",k);
		//printf("p[k]=%f\n",p[k]);
		if (a[k]!=-1){
			l_set[n_set-1]=a[k];
			res-=aw[k];
		}
		else{
			n_set--;
		}
		a[k]=u;
		aw[k]=wmax;
		res+=wmax;
	}

	g->edges=malloc(g->n1*sizeof(edge));
	j=0;
	for (i=0;i<g->n2;i++){
		if (a[i]!=-1){
			g->edges[j].u=a[i];
			g->edges[j].v=i;
			g->edges[j++].w=aw[i];
		}
	}
	g->edges=realloc(g->edges,j*sizeof(edge));
	g->e=j;
	return res;
}

void printres(bipgraph *g,char *output){
	unsigned i;
	edge ed;
	FILE* file=fopen(output,"w");
	if (g->b){
		for (i=0;i<g->e;i++) {
			ed=g->edges[i];
			fprintf(file,"%u %u %u\n",ed.v,ed.u,ed.w);
		}
	}
	else{
		for (i=0;i<g->e;i++) {
			ed=g->edges[i];
			fprintf(file,"%u %u %u\n",ed.u,ed.v,ed.w);
		}
	}
	fclose(file);
}


int main(int argc,char** argv){
	bipgraph* g;
	unsigned w;
	time_t t0,t1,t2;
	double eps=atof(argv[3]);
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);
	printf("Number of left nodes: %u\n",g->n1);
	printf("Number of right nodes: %u\n",g->n2);
	printf("Number of edges: %u\n",g->e);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Building datastructure\n");

	mkbg(g);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing optimal assignment\n");

	w=bidding(g,eps);

	printf("optimal assignment value = %u\n",w);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Writting the results in file %s\n",argv[2]);
	printres(g,argv[2]);
	freebg(g);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
