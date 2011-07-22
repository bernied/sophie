//------------------------------------------------------------------------------
/*http://www.facebook.com/careers/puzzles.php?puzzle_id=11
This program attempts to solve a problem posed by The Facebook puzzle master[1].
The problem can be summarized as follows:

You come home after a long day coding looking for some relaxation with your cat,
Sophie. She can hide in any one of multiple places in your house. It takes some
amount of time to move from hiding place to hiding place, and there is an
associated probability that she will be hiding in each particular spot. Your
goal is to find the minimum expected time to find the cat[2].

The Puzzle master gives you a list of hiding place names and the probability
that Sophie will be hiding there. You are also given a second list of hiding
place name pairs and the number of seconds it takes to move between the two
hiding places.

This specification defines a weighted graph. Each hiding place represents a node
or vertex in the graph with an associated probability (of finding the cat
there). Each pair of hiding places specifies a link from the first to the
second, and the associated time it takes to walk from the first hiding place to
the second. Note, that the graph is symmetric, so the time it takes to walk in
one direction and find the cat is equal to the time it takes to walk in the
opposite direction and find the cat. Also note the time it takes to search for
the cat at any one particular hiding place is zero (this helps keep the problem
symmetric).

In order to find Sophie we must move from hiding place to hiding place along one
of the possible links between each hiding place. This is equivalent to
traversing the nodes of the graph in some particular order. Each time we move
from hiding place to hiding place we must expend some amount of time, which is
defined by the weight of the link that was chosen to make the move.

Thus, every traversal of the graph is defined by a list of nodes and the
associated links that were used to move between the nodes. In order to
calculate the expected time for a particular traversal we must do a weighted sum
of the probability of finding Sophie at a particular node multiplied by the
total amount of time it took to traverse to this node.

The formula looks like this:

[XXX]

In the given example the Puzzle Master suggests a traversal of
front_door -> in_cabinet -> front_door -> under_bed -> behind_blinds.
We can calculate the expected time by summing the total time it takes to cross
each link times the probability that Sophie is there. In this case we get
  (0.2 front_door) * 0 seconds
+ (0.3 in_cabinet) * (0 + 2) seconds
+ (0.0 front_door) * (0 + 2 + 2) seconds
+ (0.4 under_bed) * (0 + 2 + 2 + 5) seconds
+ (0.1 behind_blinds) * (0 + 2 + 2 + 5 + 9) seconds
=
0 seconds + 0.6 seconds + 0 seconds + 3.6 seconds + 1.8 seconds
= 6 total seconds
(which just so happens to be the minimum expected time for this graph)

Three things to note here:
1. Since it takes no time to search for Sophie once we are at a node, we
can assume our starting node will require no time, hence the 0 seconds for
front_door.
2. In going from in_cabinet to front_door we are returning to a node where we
had previously checked for the cat. The probabililty for finding Sophie at a
node that we had previously been to is zero, otherwise the search would have
stopped already.
3. Crossing each link is costing us time regardless of whether that node had
been previously traversed. A previously traversed node will not add to the
nodes contribution to the expected time, but it will add to the next "new"
node's contribution to the expected time. In this case we see going from
front_door to under_bed now has a total of 9 seconds, because of the extra
traversal through front_door.

In order to simplify the problem we can make the formula link oriented instead
of node oriented. Currently for each node traversed we add its weighted
contribution to the total time. If instead we add each links contribution, we
get a formula like this:

[XXX]

For the given example we would have:
  0 seconds * (0.2 front_door + 0.3 in_cabinet + 0.0 front_door + 0.4 under_bed + 0.1 behind_blinds)
  2 seconds * (0.3 in_cabinet + 0.0 front_door + 0.4 under_bed + 0.1 behind_blinds)
+ 2 seconds * (0.0 front_door + 0.4 under_bed + 0.1 behind_blinds)
+ 5 seconds * (0.4 under_bed + 0.1 behind_blinds)
+ 9 seconds * (0.1 behind_blinds)

We can simplify the formula by eliminating all nodes that have been traversed
before. In this case we get the following for the example:

  2 seconds * (0.3 in_cabinet + 0.4 under_bed + 0.1 behind_blinds)
+ 2 seconds * (0.4 under_bed + 0.1 behind_blinds)
+ 5 seconds * (0.4 under_bed + 0.1 behind_blinds)
+ 9 seconds * (0.1 behind_blinds)

This simplifies to:

  2 seconds * (0.3 in_cabinet + 0.4 under_bed + 0.1 behind_blinds)
+(2 seconds + 5 seconds) * (0.4 under_bed + 0.1 behind_blinds)
+ 9 seconds * (0.1 behind_blinds)

Because we are dealing with probabilities we know that they will always add up
to one. Therefor we can invert each probability sum by subtracting the missing
probabilities from 1.0. If we do this we get:

[x]
  2 seconds * (1.0 - 0.2 front_door)
+(2 seconds + 5 seconds) * (1 - (0.2 front_door + 0.3 in_cabinet))
+ 9 seconds * (1 -(0.2 front_door + 0.3 in_cabinet + 0.4 under_bed))

This formula will be used as the basis for solving the problem because it
exposes the underlying structure of the problem. In order to see the structure
we must first understand the concept of a {new node ordering}.

A traversal specifies a series of nodes and a series of links that were used to
move from one node to the next. In traversing from one node to another, you can
always ask your self the question "Have I seen this node before?" In otherwords,
is this a new node or an old node. If we take a specific traversal, we can see
that it will imply a new node ordering. We can construct a list
of new nodes from the list of nodes in the traversal by adding every new node we
see in the traversal to the new node list making sure to keep the ordering the
same. 

More generally, every traversal will imply a particular new node ordering.
Multiple different traversals may imply the same new node ordering. A new node
ordering will always be no larger then the number of nodes in a particular
graph, and will be exactly that number if the traversal touches all nodes.

In the given example the traversal will imply the following new node ordering:
front_door, in_cabinet, under_bed, behind_blinds.  Another traversal that gives
the same new node ordering would be the traversal of front_door -> in_cabinet
-> front_door -> in_cabinet -> front_door -> under_bed -> behind_blinds. This
would still imply the same new node ordering.

If we look at figure [x] we can see the new node ordering appearing in each
probability term. It grows until the last summation term exposes the whole new
node ordering (excepting the last element). (Note that this formulation exposes
the starting node but hides the final node, whereas the previous formulation
exposes the last node but hides the final node).

Another advantage of this formulation is that the time term is the sum of the
seconds it takes to traverse from the previous node to the current node. This is
much more in line with how we think about graph traversals.

Lets review the steps we use to calculate the expected time for a traversal:
1. Extract the new node ordering from the traversal.
2. For each new node (except the last)
	-calculate the time it will take to go from current node to the next node
	-calculate 1.0 minus the sum of the probabilities of all the nodes traversed so far
	-multiply these two values
3. Sum the values for each node.

The general formula is as follows:
XXX

The Puzzle master is requesting the minimum expected time, not just any expected
time. Recall that every traversal implies one new node ordering, therefor there
is a many to one relationship between traversals and new node orderings. Thus
for each new node ordering there is one traversal that takes less time then all
others, ie a minimum expected time traversal[x]. This however is not gauranteed
to be the minimum expected time traversal because there are multiple new node
orderings.

We can restate the problem as follows:
  For each new node ordering find the minimum expected time traversal. Among all
  these new node orderings chose the one with the lowest expected time of all
  minimum expected time traversals. This will be the minimum expected time for
  finding Sophie.

There are two issues that need to be resolved in order to solve the problem in
this way. First we need to find all the possible new node orderings for a given
graph. Assuming the graph is complete, and that the graph contains N nodes. Then
given a fixed starting node, we can say there are (N-1)! possible new node
orderings. 

<Proof>

The second problem is to find the minimum expected time traversal for a given
new node ordering. Given a node, we know what the next new node should be. In 
order to minimize this value we must use the minimum traversal time from current
node to the destination node. This traversal may occur through one link or many.
We only care that the series of links is the minimum.

<Proof>

If a given minimum path from the current node to destination node only passes
through old nodes, there are no issues. If the minimum path does touch a new
node, then the given new node ordering is not a possible candidate for the
minimum over all expected paths.

<Proof>

Calculating all the different possible new node orderings is fairly
straightforward. However, you will end up with (N-1)! new node orderings. This
means brute force techniques are going to take too long after a fixed size.

This raises an important point. This problem of finding the minimum expected
time to find Sophie is a variation of the Traveling Salesman Problem (TSP).
The TSP attempts to find the minimum time or distance that would be required to
visit a set of cities with out visiting any city more then once. This problem
reframes the TSP in two ways.

First we say that we can visit any node more then once, because we want to find
the minimum overall time, not the minimum time if we just visited each room once.
It would not make sense to visit each room only once when we can be faster by
moving through a room that had been previously visited. 

Second the problem includes a weighting at each node. Most TSP problems talk
about weighting of links but not nodes. Having the weights at the nodes however
does not change the fundamental problem. Imagine if each node had the same
weight, then we are back to the original TSP. All the weights at the nodes will
do is alter the contribution that each node will make to the final total
expected minimum time. You can see that in the formula [x]. The early nodes get
to contribute more of their previous link cost to the final total. However,
this does not change the fact that there are still (N-1)! combinations to
search. We are still up against the wall.

Calculating the minimum expected time traversal for a particular new node
ordering will require calculating the minimum path between any given two nodes.
This can be done in O(N^3) time and O(N^2) space by using the Floyd-Warshall
algorithm. It is interesting to note that the way to solve the TSP with out the
constraint of having to visit each city once is exactly this technique of using
the minimum cost paths between each of the nodes. [See book]

In order to solve this problem I use brute force enumeration techniques up to a
given graph size (say 10 nodes). I create each possible new node ordering and
calculate its minimum expected time using formula [x]. The lowest value for each
of the new node orderings is the answer.

If there are more nodes in the graph then my threshold then the code will switch
to using a heuristic called simulated annealing. Simulated annealing will not
gaurantee the exact minimum expected time, however it can get fairly close in a
reasonable amount of time.

What begins below is the code to solve this problem. Each chunk of code will be
explained with a final explanation of how it all fits together.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>

#include "parse_cl.h"

#define MAX_BUFF 1024
#define MAX_WEIGHT UINT_MAX
#define MAX_NUMBER_VERTS 10
#define VERSION "1.0"

/*
4
front_door    .2
in_cabinet    .3
under_bed     .4
behind_blinds .1
5
front_door under_bed     5
under_bed  behind_blinds 9
front_door behind_blinds 5
front_door in_cabinet    2
in_cabinet behind_blinds 6
*/

typedef struct
{
	char* name;
	unsigned int id;
	unsigned long hash;
	float weight;
} vertex;

typedef struct
{
	vertex* source;
	vertex* dest;
	unsigned int weight;
} edge;

typedef struct
{
	unsigned int num_nodes;
	vertex* nodes;
	unsigned int num_edges;
	edge* edges;
} graph_spec;

typedef struct
{
	unsigned int num_verts;
	float* vert_weights;
	unsigned int num_edges;
	unsigned int* parent; // parent graph for path recovery
	unsigned int* edge_weights;
} graph;

vertex* get_vertex(graph_spec* spec, unsigned int id);
void print_graph(graph* g, graph_spec* spec);
void print_parent_graph(graph* g);

/*****************************************************************************/
/* Utility functions                                                         */
/*****************************************************************************/

/*
Daniel Bernstein's hash function.
See http://fr.wikipedia.org/wiki/Table_de_hachage#Fonction_de_Hachage.
*/
unsigned long
hash(char *str)
{
	unsigned long hash = 5381;
	while(*str!='\0') {
	        int c = *str;
                /* hash = hash*33 + c */
                hash = ((hash << 5) + hash) + c;
                str++;
        }
	return hash;
}

/*
Since the graph is symmetric we don't need to store the weights for the links
for both directions. We can store just the link in one direction and reverse
the source and destination in the right condition.
We borrow a trick from matrix calculations and store the array data in a format
called lower triangular. For example, a symmetric graph with five nodes (A, B,
C, D) might look like this:

  A B C D
A * 1 2 3
B 1 * 4 5
C 2 3 * 6
D 4 5 6 *

Note that along the diagonal the values are not stored. This just implies that
links from a node to itself are not relevant. When we map this graph into a
lower triangular format we get the following:

  A B C D
A * * * *
B 1 * * *
C 2 3 * *
D 4 5 6 *

The space savings comes from mapping the lower triangular array into a linear
array as such:

  1 2 3 4 5 6

In order to access a weighted link between two nodes we must map two nodes
into a specific position into the linear array. In order to do this we use
the formula
  [source * (source - 1) / 2 + destination]

*/
unsigned int
LT(unsigned int src, unsigned int dst)
{
	if (src < dst) {
		unsigned int tmp = src;
		src = dst;
		dst = tmp;
	}
	return (src * (src-1) / 2) + dst;
}

/*
find_str searches for a space delimted string in the character buffer buff. It
finds the start and ending offsets (space not included) of the string. This is
used to parse the input file.
*/
void
find_str(char* buff, unsigned int max, unsigned int offset, unsigned int* start, unsigned int* end)
{
	unsigned int i = offset;
	while(buff[i] == ' ' || buff[i] == '\t') {
		i++;
	}
	*start = i;
	while(i < max && (buff[i] != ' ' && buff[i] != '\t')) {
		i++;
	}
	*end = i;
}

char*
substring(char* buff, unsigned int start, unsigned int end)
{
	char* str;
	unsigned int len = end - start;
	str = (char*) malloc(len+1);
	if (!str) {
		return 0;
	}
	strncpy(str, buff+start, len);
	str[len] = '\0';

	return str;
}

/* taken from http://www.dreamincode.net/code/snippet4206.htm */
/*
To generate a permuation we use the following algorithm...
*/
int
lex_permute(unsigned int *array, int size) {
  int l, r, tmp;
  if (size < 2)
    return 0;
  l = size - 2;
  while (l >= 0 && array[l] > array[l+1])
    --l;

  if (l < 0)
    return 0;
  r = size - 1;
  while (array[r] < array[l])
    --r;
  tmp = array[l];
  array[l] = array[r];
  array[r] = tmp;
  r = size - 1;
  ++l;
  while (l < r) {
    tmp = array[l];
    array[l++] = array[r];
    array[r--] = tmp;
  }
  return 1;
}

unsigned int random_int(unsigned int min, unsigned int max)
{
	return rand() % (max - min + 1) + min;
}

double random_float()
{
	return rand() / (((double)RAND_MAX)+1.0);
}
/*****************************************************************************/
/* graph - An ADT to contain the minimum distance between all pairs          */
/*****************************************************************************/

/* Floyd-Warshall algorithm from "The Algorithm Design Manual", p. 211, w/ minor tweaks. */
void
all_pairs_shortest_paths_with_parents(graph * g)
{
	unsigned int i, j;
	unsigned int k;
	unsigned long through_k;
	register unsigned int nv = g->num_verts;

	for (i=0; i < nv; i++) {
		for (j=0; j < nv; j++) {
			if (i == j) continue;
//			if (g->edge_weights[LT(i, j)] != MAX_WEIGHT)
			{
				g->parent[(i * nv) + j] = i;
			}
		}
	}

	// This could be optimized to do 1/2 the work! How?
	for (k=0; k < nv; k++) {
		for (i=0; i < nv; i++) {
			for (j=0; j < nv; j++) {
				if (i == j) continue; // if equal then link to self!
				through_k = g->edge_weights[LT(i,k)] + g->edge_weights[LT(k,j)];
				if (through_k <= g->edge_weights[LT(i,j)]) { // use <= to make sure parent array is set!
					g->edge_weights[LT(i,j)] = through_k;
					g->parent[(i * nv) + j] = g->parent[(k * nv) + j];
				}
			}
		}
	}
}

void
all_pairs_shortest_paths(graph * g)
{
	unsigned int i, j;
	unsigned int k;
	unsigned long through_k;
	register unsigned int nv = g->num_verts;

	// This could be optimized to do 1/2 the work! How?
	for (k=0; k < nv; k++) {
		for (i=0; i < nv; i++) {
			for (j=0; j < nv; j++) {
				if (i == j) continue; // if equal then link to self!
				through_k = g->edge_weights[LT(i,k)] + g->edge_weights[LT(k,j)];
				if (through_k <= g->edge_weights[LT(i,j)]) { // use <= to make sure parent array is set!
					g->edge_weights[LT(i,j)] = through_k;
				}
			}
		}
	}
}

/*
http://www.cs.cornell.edu/~wdtseng/icpc/notes/graph_part3.pdf
The meaning of parent[i][j] is this:
"For some shortest path from i to j, the vertex right before j is parent[i][j]."
We initialize all entries parent[i][j] to i when we have an edge, and if not we set it to 1.
The update process of the parent array basically mean the following:
"If going from i to j through k (i -> k -> j) is an improvement, then we set the parent of j in
the new shortest path to the parent of the path k -> j."
To recover the path from i to j, simply recurse on parent[i][j]'s entries.
*/
int
extract_min_path(graph* g, unsigned int src, unsigned int dst, unsigned int* path, unsigned int* size)
{
	if (path == NULL || *size < 2 || src == dst) {
		return 0;
	}
	
	unsigned int i=0;
	do {
		path[i++] = src;
		src = g->parent[(g->num_verts * src) + dst];
		if (src == MAX_WEIGHT) {
			printf("unable to extract min path!\n");
			return 0;
		}
	} while(src != path[i-1] && i < *size);
	path[i] = dst;

	*size = i+1;

	return 1;
}

graph*
init_graph(graph_spec* spec)
{
	unsigned int i;

	if (!spec) {
		return NULL;
	}

	graph* g = (graph*) malloc(sizeof(graph));
	if (!g) {
		return NULL;
	}

	unsigned int count =  spec->num_nodes * (spec->num_nodes - 1) / 2;
	g->edge_weights = (unsigned int*) malloc(count * sizeof(unsigned int));
	if (!g->edge_weights) {
		free(g);
		return NULL;
	}
	for (i=0; i < count; i++) {
		g->edge_weights[i] = MAX_WEIGHT;
	}

	unsigned long parent_size = spec->num_nodes * spec->num_nodes;
	g->parent = (unsigned int*) malloc(parent_size * sizeof(unsigned int)); // LAMb: could make this a sparse graph
	if (!g->parent) {
		free(g->edge_weights);
		free(g);
		return NULL;
	}
	for (i=0; i < parent_size; i++) {
		g->parent[i] = MAX_WEIGHT;
	}

	g->vert_weights = (float*) malloc(spec->num_nodes * sizeof(float));
	if (!g->vert_weights) {
		free(g->edge_weights);
		free(g->parent);
		free(g);
		return NULL;
	}

	g->num_verts = spec->num_nodes;
	g->num_edges = spec->num_edges;

	for(i=0; i < spec->num_nodes; i++) {
		g->vert_weights[spec->nodes[i].id] = spec->nodes[i].weight;
	}

	for(i=0; i < spec->num_edges; i++) {
		edge* e = &spec->edges[i];
		g->edge_weights[LT(e->source->id, e->dest->id)] = e->weight;
	}

	all_pairs_shortest_paths(g);

	//print_graph(g, spec);
	//print_parent_graph(g);

	return g;
}

int
is_graph_disconnected(graph* g)
{
	if (!g) {
		return 0;
	}
	
	for (unsigned int i=0; i < g->num_edges; i++) {
		switch(g->edge_weights[i]) {
			case 0:
			case MAX_WEIGHT:
				return 1;
		}
	}
	
	return 0;
}

void
free_graph(graph* g)
{
	if (g) {
		if (g->vert_weights) {
			free(g->vert_weights);
		}
		if (g->edge_weights) {
			free(g->edge_weights);
		}
		if (g->parent) {
			free(g->parent);
		}
		free(g);
	}
}

void
print_parent_graph(graph* g)
{
	unsigned int i, j;

	for (i=0; i < g->num_verts; i++) {
		for (j=0; j < g->num_verts; j++) {
			printf("%u,\t", g->parent[(i * g->num_verts) + j]);
		}
		printf("\n");
	}
}

void
print_graph(graph* g, graph_spec* spec)
{
	unsigned int i, j;
	vertex *src, *dst;

	if (spec) {
		for(i=0; i < g->num_verts; i++) {
			for (j=0; j < g->num_verts; j++) {
				if (i != j) {
					src = get_vertex(spec, i);
					dst = get_vertex(spec, j);
					if (src != NULL && dst != NULL) {
						printf("%s -> %s = %u\n", src->name, dst->name, g->edge_weights[LT(i, j)]);
					}
				}
			}
		}
	}
	else {
		for(i=0; i < g->num_verts; i++) {
			for (j=0; j < g->num_verts; j++) {
				if (i != j) {
					printf("%u -> %u = %u\n", i, j, g->edge_weights[LT(i, j)]);
				}
			}
		}
	}
}

void
print_tour(graph_spec* spec, unsigned int* tour, unsigned int size, float time)
{
	vertex* v;
	for (unsigned int i=0; i < size; i++) {
		v = get_vertex(spec, tour[i]);
		if (v) {
			printf("%s -> ", v->name);
		}
		else {
			printf("%u -> ", tour[i]);
		}
	}
	printf("%0.2f\n", time); // LAMb: print 2 values after decimal point!
}

float
calculate_probability(graph* g, unsigned int* tour, unsigned int pos)
{
	return g->vert_weights[tour[pos]];
}

float
calculate_time(graph* g, unsigned int* tour, unsigned int tour_size)
{
	if (tour == NULL) {
		return -1.0f;
	}
	float probability = 1.0f;
	float time = 0.0;
	unsigned int weight;

	for (unsigned int i=0; i < tour_size-1; i++) {
		probability -= calculate_probability(g, tour, i);
		weight = g->edge_weights[LT(tour[i], tour[i+1])];
		time += ((float)weight) * probability;
	}

	return time;
}

int
is_valid_tour(graph* g, graph_spec* gs, unsigned int* tour, unsigned int size)
{
	if (tour == NULL) {
		return 0;
	}
	unsigned int path_size, i, j;
	char* node_visited = (char*) calloc(size, sizeof(char));
	if (!node_visited) {
		return -1;
	}
	unsigned int* min_path = (unsigned int*) malloc(size * sizeof(unsigned int));
	if (!min_path) {
		return -1;
	}

	node_visited[tour[0]] = 1;
	for (i=0; i < size-1; i++) {
		path_size = size;
		if (!extract_min_path(g, tour[i], tour[i+1], min_path, &path_size)) {
			return -1;
		}
		for (j=0; j < path_size-1; j++) {
			if (!node_visited[min_path[j]]) {
//				print_tour(gs, min_path, path_size, 0.0f);
				return 0;
			}
		}
		if (node_visited[min_path[j]]) {
			return -1;
		}
		node_visited[min_path[j]] = 1; // visit last node because it should be new!
	}
	return 1;
}

float
brute_force_calc(graph* g, graph_spec* spec, struct arg_t* my_args)
{
	unsigned int* permutation = (unsigned int*) malloc(g->num_verts * sizeof(unsigned int));
	if (permutation == NULL) {
		return -1.0; // LAMb: lame!
	}	
	for (unsigned int i=0; i < g->num_verts; i++) {
		permutation[i] = i;
	}

	float time;
	float minimum = calculate_time(g, permutation, g->num_verts);
	if (my_args->p) {
        print_tour(spec, permutation, g->num_verts, minimum);
    }
	while(lex_permute(permutation+1, g->num_verts-1)) {
		/* if (!is_valid_tour(g, spec, permutation, g->num_verts)) { // LAMb: returns -1! */
		/* 	continue; */
		/* } */
		time = calculate_time(g, permutation, g->num_verts);
		if (time < minimum) {
            if (my_args->p) {
                print_tour(spec, permutation, g->num_verts, time);
            }
			minimum = time;
		}
	}

	free(permutation);
	return minimum;
}

void
swap(unsigned int* p, unsigned int nv, unsigned int i1, unsigned int i2)
{
	unsigned int temp = p[i1];
	p[i1] = p[i2];
	p[i2] = temp;
}

int
transition(graph* g, unsigned int* permutation, unsigned int nv, unsigned int i1, unsigned int i2)
{
	if (i1 == i2) {
		return 0;
	}

	swap(permutation, nv, i1, i2);
//	if (!is_valid_tour(g, NULL, permutation, nv)) {
//		swap(permutation, nv, i1, i2);
//		return 0;
//		}

	return 1;
}

/* Taken from "The Algorithm Design Manual", p. 257-258, w/ minor tweaks. */
#define E 2.718
unsigned int
anneal(graph* g, unsigned int* permutation, struct arg_t* my_args)
{
	int i1, i2;
	int i, j;
	double temperature;
	double current_cost, start_cost, cost;
	double delta;
	double merit, flip;
	double exponent;
	unsigned int solutions =0;
	unsigned int nv = g->num_verts;

	temperature = my_args->i;
	current_cost = calculate_time(g, permutation, nv);

	for (i=1; i <= my_args->s; i++) {
		temperature *= my_args->f;
		start_cost = current_cost;

		for(j=1; j <= my_args->t; j++) {
			//prev_cost = calculate_time(g, permutation, nv);
			do {
				i1 = random_int(1, nv-1); // leave 1st position in tact!
				i2 = random_int(1, nv-1);
			} while(!transition(g, permutation, nv, i1, i2));
			cost = calculate_time(g, permutation, nv);
			delta = cost - current_cost;

//printf("delta = %f current_cost = %f cost = %f temp = %f\n", delta, current_cost, cost, temperature);

			if (delta < 0) {
				current_cost = cost;
			}
			else {
				flip = random_float();
				exponent = (-delta/current_cost)/(my_args->k*temperature);
				merit = exp(exponent);
//printf("merit = %f  flip=%f  exponent=%f\n",merit,flip,exponent);

				if (merit > flip) {
					current_cost = current_cost + delta;
				}
				else {
					swap(permutation, nv, i1, i2);
				}
			}
			solutions++;
		}

		if ((current_cost - start_cost) < 0.0) {
			temperature = temperature / my_args->f;
		}
	}
	return solutions;
}

float
heuristic_calc(graph* g, graph_spec* spec, struct arg_t* my_args)
{
	unsigned int* permutation = (unsigned int*) malloc(g->num_verts * sizeof(unsigned int));
	unsigned int* best_perm = NULL;
	if (permutation == NULL) {
		return -1.0; // LAMb: lame!
	}	
	if (my_args->e) {
		best_perm = (unsigned int*) malloc(g->num_verts * sizeof(unsigned int));
	}
	unsigned int nv = g->num_verts;
	unsigned int i1, i2;

	float total_time = FLT_MAX;
	for (int i=0; i < my_args->r; i++) {
		if (i > 0 && my_args->e && best_perm) {
			for (unsigned int j=0; j < nv; j++) {
				permutation[j] = best_perm[j];
			}
		}
		else {
			for (unsigned int j=0; j < nv; j++) {
				permutation[j] = j;
			}
		}

		if (my_args->z) {
			for (unsigned int k=0; k < nv; k++) {
				i1 = random_int(1, nv-1); // leave 1st position in tact!
				i2 = random_int(1, nv-1);
				transition(g, permutation, nv, i1, i2);
			}
		}
		
		anneal(g, permutation, my_args);
		
		float this_time = calculate_time(g, permutation, nv);
		if (this_time < total_time) {
			total_time = this_time;
			if (my_args->e && best_perm) {
				for (unsigned int j=0; j < nv; j++) {
					best_perm[j] = permutation[j]; // LAMb: use copy function?
				}
			}
		}
	}
//	print_tour(spec, permutation, g->num_verts, total_time);
	free(permutation);
	if (best_perm)
		free(best_perm);

	return total_time;
}

/*****************************************************************************/
/* graph_spec - ADT that mirrors the specification of the problem file.      */
/*****************************************************************************/

vertex*
get_vertex(graph_spec* spec, unsigned int id)
{
	for(unsigned int i=0; i < spec->num_nodes; i++) {
		if (spec->nodes[i].id == id) {
			return &spec->nodes[i];
		}
	}
	return NULL;
}

vertex*
find_vertex(vertex* v, unsigned int num_verticies, char* name)
{
	if (!v || !name) {
		return NULL;
	}
	unsigned long h = hash(name);

	for (unsigned int i=0; i < num_verticies; i++) {
		if (v[i].hash == h) {
			if (strcmp(v[i].name, name) == 0) {
				return v + i;
			}
		}
	}
	return NULL;
}

int
init_vertex(char* line, unsigned int num_chars, unsigned int id, vertex* v)
{
	if (!v) {
		return 0;
	}

	unsigned int start;
	char* name;
	char* weight_str;
	float weight;

	unsigned int i=0;
	find_str(line, num_chars, i, &start, &i);
	name = substring(line, start, i);
	if (!name) {
		return 0;
	}

	find_str(line, num_chars, i, &start, &i);
	weight_str = substring(line, start, i);
	if (!weight_str) {
		free(name);
		return 0;
	}
	weight = (float) atof(weight_str);
	free(weight_str);
	
	v->id = id;
	v->hash = hash(name);
	v->name = name;
	v->weight = weight;

	return 1;
}

int
init_edge(char* line, unsigned int num_chars, vertex* verticies, unsigned int num_verticies, edge* e)
{
	if (!e) {
		return 0;
	}

	unsigned int start;
	char* name;
	char* weight_str;
	unsigned int weight;

	unsigned int i=0;
	find_str(line, num_chars, i, &start, &i);
	name = substring(line, start, i);
	if (!name) {
		return 0;
	}
	e->source = find_vertex(verticies, num_verticies, name);
	free(name);
	if (!e->source) {
		return 0;
	}

	find_str(line, num_chars, i, &start, &i);
	name = substring(line, start, i);
	if (!name) {
		return 0;
	}
	e->dest = find_vertex(verticies, num_verticies, name);
	free(name);
	if (!e->dest) {
		return 0;
	}

	find_str(line, num_chars, i, &start, &i);
	weight_str = substring(line, start, i);
	if (!weight_str) {
		return 0;
	}

	weight = atoi(weight_str);
	free(weight_str);
	e->weight = weight;

	return 1;
}

void
free_verticies(vertex* v, unsigned int num_verticies)
{
	for (unsigned int i=0; i < num_verticies; i++) {
		free(v[i].name);
	}
	free(v);
}

void
free_edges(edge* e, unsigned int num_edges)
{
	// Assumes free_verticies is called to free vertex pointers.
	free(e);
}

void
free_graph_spec(graph_spec* g)
{
	if (!g) {
		return;
	}
	if (g->nodes) {
		free_verticies(g->nodes, g->num_nodes);
	}
	if (g->edges) {
		free(g->edges);
	}
	free(g);
}

/*****************************************************************************/
/* File I/O - Read the input file to create a graph_spec.                    */
/*****************************************************************************/

unsigned int
read_line(FILE* f, char* buff, unsigned int max)
{
	fpos_t position;
	unsigned int i =0;
	int c;

	do {
		c = fgetc(f);
		switch(c) {
			case '\r':
				fgetpos(f, &position);
				c = fgetc(f);
				if (c != '\n') {
					fsetpos(f, &position);
				}
			case '\n':
			case -1:
				return i;
			default:
				if (i < max) {
					buff[i++] = (char) c;
				}
				else {
					return i;
				}
		}
	} while(1);
}

vertex*
read_verticies(FILE* f, unsigned int* num)
{
	char line[MAX_BUFF];
	
	unsigned int num_chars = read_line(f, line, MAX_BUFF);
	line[num_chars] = 0;
	unsigned int n = (unsigned int) atol(line);
	
	vertex* verticies = (vertex*) malloc(sizeof(vertex) *  n);

	for(unsigned int i=0; i < n; i++) {
		num_chars = read_line(f, line, MAX_BUFF);
		if (!init_vertex(line, num_chars, i, verticies + i)) {
			free_verticies(verticies, i);
			return NULL;
		}
	}

	*num = n;
	return verticies;
}

edge*
read_edges(FILE* f, vertex* verticies, unsigned int num_verticies, unsigned int* num)
{
	char line[MAX_BUFF];
	
	unsigned int num_chars = read_line(f, line, MAX_BUFF);
	line[num_chars] = 0;
	unsigned int n = (unsigned int) atol(line);
	
	edge* edges= (edge*) malloc(sizeof(edge) *  n);

	for(unsigned int i=0; i < n; i++) {
		num_chars = read_line(f, line, MAX_BUFF);
		if (!init_edge(line, num_chars, verticies, num_verticies, edges + i)) {
			free_edges(edges, i);
			return NULL;
		}
	}

	*num = n;
	return edges;
}

graph_spec*
read_graph_spec(FILE* file)
{
	graph_spec* g = (graph_spec*) calloc(1, sizeof(graph_spec));
	if (!g) {
		return NULL;
	}
	
	g->nodes = 	read_verticies(file, &g->num_nodes);
	if (!g->nodes) {
		free(g);
		return NULL;
	}

	g->edges = read_edges(file, g->nodes, g->num_nodes, &g->num_edges);
	if (!g->edges) {
		free_graph_spec(g);
		return NULL;
	}

	return g;
}

#define INITIAL_TEMP 1.0
#define COOLING_STEPS 500
#define COOLING_FRACTION 0.97
#define STEPS_PER_TEMP 1000
#define NUM_RUNS 1
#define K_FACTOR 0.01
void
init_default_args(struct arg_t* my_args)
{
    my_args->h = 0;
    my_args->V = 0;
    my_args->p = 0;
    my_args->i = INITIAL_TEMP;
    my_args->t = STEPS_PER_TEMP;
    my_args->s = COOLING_STEPS;
    my_args->f = COOLING_FRACTION;
    my_args->k = K_FACTOR;
    my_args->r = NUM_RUNS;
    my_args->v = MAX_NUMBER_VERTS;
    my_args->optind = 0;
}

int
main(int argc, char* argv[])
{
	//if (argc != 2) {
	//	printf("Expecting single file name as input.");
	//	exit(EXIT_FAILURE);
	//}
	float total_time;
	FILE * f;
    struct arg_t my_args;
    
    init_default_args(&my_args);
    Cmdline(&my_args, argc, argv);
    if (argc == 1 || my_args.optind+1 != argc) {
        usage(-1, argv[0]);
    }
    
    if (my_args.V) {
        printf("sophie version %s\n", VERSION);
    }

	f = fopen((const char*) argv[my_args.optind], "r");
	if (!f) {
		printf("Unable to open file %s.", argv[1]);
		exit(EXIT_FAILURE);
	}
	graph_spec* gs = read_graph_spec(f);
	fclose(f);
	if (!gs) {
		printf("Input file invalid, failed to process.\n");
		exit(EXIT_FAILURE);
	}
	graph* g = init_graph(gs);
	if (!g) {
		printf("Failed to initialize the graph.\n");
		exit(EXIT_FAILURE);
	}
	if (is_graph_disconnected(g)) {
		printf("-1.00\n");
		exit(EXIT_SUCCESS);
	}

	srand(time(0));
	if (my_args.b || g->num_verts <= my_args.v) {
		total_time = brute_force_calc(g, gs, &my_args);
	}
	else {
		total_time = heuristic_calc(g, gs, &my_args);
	}

	free_graph_spec(gs);
	free_graph(g);

	printf("%0.2f\n", total_time); // LAMb: print 2 values after decimal point!
	exit(EXIT_SUCCESS);
}
