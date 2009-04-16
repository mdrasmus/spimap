#include <math.h>
#include "Tree.h"
#include "phylogeny.h"
#include "birthdeath.h"


namespace spidir
{

extern "C" {


// returns the number of labeled histories with 'ngenes' surviving lineages
int inumHistories(int ngenes)
{
    // gaurd against overflow
    assert(ngenes <= 9);

    int n = 1;
    for (int i=2; i<=ngenes; i++) {
        n *= i*(i-1) / 2;
    }
    return n;
}


// returns the number of labeled histories with 'ngenes' surviving lineages
double numHistories(int ngenes)
{
    double n = 1;
    for (int i=2; i<=ngenes; i++) {
        n *= i*(i-1) / 2;
    }
    return n;
}

// returns the number of labeled histories exist for the given tree topology
// NOTE: assumes binary tree
int inumTopologyHistories(Tree *tree)
{
    int n = 1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];

        if (node->isLeaf()) {
            ninternals[node->name] = 0;
        } else {
            // count internal children
            const int right = ninternals[node->children[0]->name];
            const int left = ninternals[node->children[1]->name];
            ninternals[node->name] = 1 + right + left;
            n *= choose(right + left, right);
        }
    }

    return n;
}


// returns the number of labeled histories exist for the given tree topology
// NOTE: assumes binary tree
double numTopologyHistories(Tree *tree)
{
    double n = 1;

    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];

        if (node->isLeaf()) {
            ninternals[node->name] = 0;
        } else {
            // count internal children
            const int right = ninternals[node->children[0]->name];
            const int left = ninternals[node->children[1]->name];
            ninternals[node->name] = 1 + right + left;
            n *= fchoose(right + left, right);
        }
    }

    return n;
}


// returns the number of labeled histories exist for the given tree topology
// uses the subtree starting at root and going until leaves.
// NOTE: assumes binary tree
double numSubtopologyHistories(Tree *tree, Node *root, ExtendArray<Node*> &leaves)
{
    double n = 1;
    
    // get nodes in post order
    ExtendArray<Node*> queue(0, tree->nnodes);

    // count number of descendant internal nodes
    ExtendArray<int> visited(tree->nnodes);
    for (int i=0; i<visited.size(); i++)
        visited[i] = 0;

    // count number of descendant internal nodes
    ExtendArray<int> ninternals(tree->nnodes);

    // process leaves
    for (int i=0; i<leaves.size(); i++) {
        Node *node = leaves[i];
        ninternals[node->name] = 0;
        
        // queue parent
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    
    // go up tree until root
    for (int i=0; i<queue.size(); i++) {
        Node *node = queue[i];
        
        // do not process a node until both children are processed
        if (visited[node->name] != 2)
            continue;
        
        // count internal children
        const int right = ninternals[node->children[1]->name];
        const int left = ninternals[node->children[0]->name];
        ninternals[node->name] = 1 + right + left;
        n *= fchoose(right + left, right);

        if (node == root)
            return n;

        visited[node->name]++;
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    // Note: this will occur for subtree like
    //    root
    //     |
    //     +
    //    / \  
    //   /   \  
    //  A     B
    // 
    return n;
}

int intCmp(const void *_a, const void *_b)
{
    return *((int*) _a) - *((int*) _b);
}


int numRedunantTopologies_helper(Node *node, 
				 bool *leafset, 
				 int *hashids)
{
    if (leafset[node->name] == true)
	return 0;
    else {
	if (node->nchildren == 1)
	    return numRedunantTopologies_helper(node->children[0],
						leafset, hashids);
	else if (node->nchildren == 2) {
	    int n = int(hashids[node->children[0]->name] == 
			hashids[node->children[1]->name]);
	    return n + numRedunantTopologies_helper(node->children[0],
						    leafset, hashids)
		     + numRedunantTopologies_helper(node->children[1],
						    leafset, hashids);
	} else {
	    // NOTE: cannot handle multifurcating nodes
	    assert(0);
	}
    }
}

double numRedunantTopologies(Tree *tree, Node *root, 
			     ExtendArray<Node*> &leaves, 
			     int *hashids)
{

    // get nodes in post order
    ExtendArray<Node*> queue(0, 2 * tree->nnodes);

    // initialize visited array to zeros
    ExtendArray<int> visited(tree->nnodes);
    for (int i=0; i<visited.size(); i++)
        visited[i] = 0;
    
    // process leaves
    for (int i=0; i<leaves.size(); i++) {
        Node *node = leaves[i];        
        // queue parent
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    int nmirrors = 0;
    
    // go up tree until root
    for (int i=0; i<queue.size(); i++) {
        Node *node = queue[i];
        
        // do not process a node until all children are processed
        if (visited[node->name] != node->nchildren)
            continue;
        
        // count internal children
        if (node->nchildren == 2) {
            nmirrors += int(hashids[node->children[0]->name] == 
                            hashids[node->children[1]->name]);
        } else {
            // we do not handle multifurcating nodes
            assert(node->nchildren < 2);
        }

        if (node == root)
            break;

        visited[node->name]++;
        visited[node->parent->name]++;
        queue.append(node->parent);
    }

    //printf("queue.size = %d\n", queue.size());


    // get hashes
    ExtendArray<int> leafhashes(0, leaves.size());
    for (int i=0; i<leaves.size(); i++)
	leafhashes.append(hashids[leaves[i]->name]);
    
    qsort((void*) leafhashes.get(), leafhashes.size(), sizeof(int), intCmp);

    double val = 1.0;
    double colorsize = 1;
    for (int i=1; i<leaves.size(); i++) {
	if (leafhashes[i] != leafhashes[i-1]) {
	    // val *= factorial(colorsize)
	    for (double j=2; j<=colorsize; j+=1.0)
		val *= j;
	    colorsize = 1.0;
	} else {
	    colorsize += 1.0;
	}
    }
    for (double j=2; j<=colorsize; j+=1.0)
	val *= j;

    //printf("c val=%f, nmirrors = %d\n", val, nmirrors);
    for (int i=0; i<nmirrors; i++)
        val /=  2.0;
    //printf("val=%f\n", val);
    return val;
}



double numRedunantTopologies2(Tree *tree, Node *root, 
			     ExtendArray<Node*> &leaves, 
			     int *hashids)
{

    // get hashes
    ExtendArray<int> leafhashes(0, leaves.size());
    for (int i=0; i<leaves.size(); i++)
	leafhashes.append(hashids[leaves[i]->name]);
    
    qsort((void*) leafhashes.get(), leafhashes.size(), sizeof(int), intCmp);

    double val = 1.0;
    double colorsize = 1;
    for (int i=1; i<leaves.size(); i++) {
	if (leafhashes[i] != leafhashes[i-1]) {
	    // val *= factorial(colorsize)
	    for (double j=2; j<=colorsize; j+=1.0)
		val *= j;
	    colorsize = 1.0;
	} else {
	    colorsize += 1.0;
	}
    }
    for (double j=2; j<=colorsize; j+=1.0)
	val *= j;

    double num = 1.0;
    for (double j=2; j<=leaves.size(); j+=1.0)
	num *= j;

    return num / val;
}



float birthDeathTopology(Node *node, float birthRate, float deathRate,
                         ExtendArray<Node*> &leaves)
{
    // numTopologyHistories(subtree) / numHistories(ngenes)

    return numHistories(leaves.size());
}

// computes the entries of the doom probabilty table
void calcDoomTable(Tree *tree, float birthRate, float deathRate, int maxdoom,
                   float *doomtable)
{
    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];
        
        if (node->isLeaf()) {
            doomtable[node->name] = -INFINITY;
        } else {
            float prod = 0.0;
            
            for (int j=0; j<node->nchildren; j++) {
                Node *child = node->children[j];
                float sum = 0.0;

                for (int ndoom=0; ndoom<=maxdoom; ndoom++) {                    
                    sum += (birthDeathCount(ndoom, child->dist, 
                                            birthRate, deathRate) *
                            ipow(exp(doomtable[child->name]), ndoom));
                }

                prod += log(sum);
            }
	    
            doomtable[node->name] = prod;
        }
    }
}



int countDups(Node *node, int *events) 
{
    int count = 0;
    if (events[node->name] == EVENT_DUP) {
        count++;
    
        for (int i=0; i<node->nchildren; i++)
            count += countDups(node->children[i], events);
    }
    
    return count;
}


// stores the leaves of the subtree below node that reconciles to snode
void getSpecSubtree(Node *node, Node *snode, int *recon, int *events,
                    ExtendArray<Node*> &nodes)
{
    for (int i=0; i<node->nchildren; i++) {
        Node *child = node->children[i];

        // only consider nodes the reconcile to snode
        if (recon[child->name] == snode->name) {
            if (events[child->name] == EVENT_SPEC ||
                events[child->name] == EVENT_GENE)
            {
                nodes.append(child);
            } else {
                getSpecSubtree(child, snode, recon, events, nodes);
            }
        }
    }
}



class KeyList {
public:
    KeyList(int key) : 
        key(key),
        next(NULL)
    {}

    int key;
    KeyList *next;
};

class HashNode {
public:
    HashNode(int nodeid, KeyList *start, KeyList *end, int len) :
        nodeid(nodeid),
	start(start),
	end(end),
	len(len)
    {}

    int nodeid;
    KeyList *start;
    KeyList *end;
    int len;
};


int hashNodeCmp(const void *_a, const void *_b)
{
    HashNode *a = *((HashNode**) _a);
    HashNode *b = *((HashNode**) _b);
    
    // first compare diffs
    int diff = a->len - b->len;
    if (diff) {
    	return diff;
    } else {
        // compare keys
	KeyList *keya = a->start;
	KeyList *keyb = b->start;

	while (true) {
	    diff = keya->key - keyb->key;
	    if (diff)
		return diff;
            if (keya == a->end || keyb == b->end)
                break;
            keya = keya->next;
            keyb = keyb->next;
	}

	// both hashes are the same
	// 1. same length
	// 2. same key subsequence
	return 0;
    }
}



void getHashIds(Tree *tree, int *recon, int *hashids)
{

    ExtendArray<HashNode*> hashnodes(tree->nnodes);
    ExtendArray<KeyList*> keylist(tree->nnodes);

    // get post order of nodes
    ExtendArray<Node*> postnodes(0, tree->nnodes);
    getTreePostOrder(tree, &postnodes);    

    // build hash nodes
    for (int i=0; i<postnodes.size(); i++)
    {
        Node *node=postnodes[i];
        
        if (node->isLeaf()) {
            KeyList *key = new KeyList(recon[node->name]);
            keylist[node->name] = key;
            hashnodes[node->name] = new HashNode(node->name, key, key, 1);
        } else {
            if (node->nchildren == 1) {                
                KeyList *key = new KeyList(-1);
                keylist[node->name] = key;
                HashNode *hnode1 = hashnodes[node->children[0]->name];

                // join lists: list1 -> key
                hashnodes[node->name] = 
                    new HashNode(node->name, hnode1->start, key, 
                                 hnode1->len + 1);
                hnode1->end->next = key;
                
            } else if (node->nchildren == 2) {
                KeyList *key = new KeyList(-2);
                keylist[node->name] = key;
                HashNode *hnode1 = hashnodes[node->children[0]->name];
                HashNode *hnode2 = hashnodes[node->children[1]->name];
                int len = hnode1->len + hnode2->len + 1;
                int cmp = hashNodeCmp(&hnode1, &hnode2);

                if (cmp <= 0) {
                    // join lists: list1 -> list2 -> key
                    hashnodes[node->name] = new HashNode(node->name,
                                                         hnode1->start, 
                                                         key,
                                                         len);
                    hnode1->end->next = hnode2->start;
                    hnode2->end->next = key;
                } else {
                    // join lists: list2 -> list1 -> key
                    hashnodes[node->name] = new HashNode(node->name,
                                                         hnode2->start, 
                                                         key,
                                                         len);
                    hnode2->end->next = hnode1->start;
                    hnode1->end->next = key;
                }
            } else {
                // cannot handle multifurcating nodes
                assert(0);
            }
        }
    }

    // sort hashnodes
    qsort((void*) hashnodes.get(), hashnodes.size(), sizeof(HashNode*),
	  hashNodeCmp);
    
    int hashid = 0;
    hashids[hashnodes[0]->nodeid] = hashid;
    for (int i=1; i<hashnodes.size(); i++) {
	// use new hashid if nodes differ
	if (hashNodeCmp(&hashnodes[i], &hashnodes[i-1]))
	    hashid++;
	hashids[hashnodes[i]->nodeid] = hashid;
    }
    

    // clean up
    for (int i=0; i<tree->nnodes; i++) {
	delete hashnodes[i];
        delete keylist[i];
    }

    //printf("hashids = ");
    //printIntArray(hashids, tree->nnodes);
    //printf("\n");
}



// TODO: does not handle branches above the species tree root yet
// NOTE: assumes binary species tree
float birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birthRate, float deathRate,
                          float *doomtable, int maxdoom)
{

    double prob = 0.0;
    ExtendArray<Node*> subleaves(0, tree->nnodes);
    
    // catch undefined params
    if (birthRate == deathRate)
        deathRate = .99 * birthRate;
    
    ExtendArray<int> hashids(tree->nnodes);
    getHashIds(tree, recon, hashids);

    // preroot duplications
    //if (events[tree->root->name] == EVENT_DUP)

    // loop through speciation nodes in tree
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = tree->nodes[i];        
        if (events[node->name] == EVENT_SPEC) {

            // loop through nodes u \in child(R(v))
            Node *snode = stree->nodes[recon[node->name]];
            for (int j=0; j<snode->nchildren; j++) {
                Node *schild = snode->children[j];

                // get subtree that reconciles to snode
                subleaves.clear();
                getSpecSubtree(node, schild, recon, events, subleaves);

		double nhist, thist;

		if (subleaves.size() == 0) {
		    nhist = 1.0;
		    thist = 1.0;
		} else {
		    nhist = numSubtopologyHistories(tree, node, subleaves) *
			numRedunantTopologies2(tree, node, subleaves, hashids);
		    thist = numHistories(subleaves.size());		    

		    // correct subtrees that have leaves
		    if (subleaves[0]->isLeaf()) {
			nhist *= numRedunantTopologies(tree, node, 
						       subleaves, 
						       hashids);
		    }
		}
		
                // sum over ndoom
                double sum = 0.0;
                for (int ndoom=0;  ndoom<=maxdoom; ndoom++) {
                    sum += (birthDeathCount(subleaves.size() + ndoom, 
                                            schild->dist,
                                            birthRate, deathRate) *
                            ipow(exp(doomtable[schild->name]), ndoom));
                }

                prob += log(nhist) - log(thist) + log(sum);
                //printf("c prod2 = %f, %d\n", sum, subleaves.size());
            }
        }
    }

    //printf("c prob = %f\n", prob);

    ExtendArray<Node*> leaves(0, tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
	if (tree->nodes[i]->isLeaf())
	    leaves.append(tree->nodes[i]);
    }
    double x = numRedunantTopologies(tree, tree->root, leaves, hashids);
    prob -= log(x);
    //printf("c x = %f\n", x);
    
    return prob;
}



// Convenience function
// Adds and removes implied species nodes to the gene tree
// NOTE: assumes binary species tree
float birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birthRate, float deathRate,
                              float *doomtable, int maxdoom)
{
    ExtendArray<int> recon2(0, tree->nnodes);
    recon2.extend(recon, tree->nnodes);

    ExtendArray<int> events2(0, tree->nnodes);
    events2.extend(events, tree->nnodes);


    int addedNodes = addImpliedSpecNodes(tree, stree, recon2, events2);
    float p = birthDeathTreePrior(tree, stree, recon2, events2, 
                                  birthRate, deathRate,
                                  doomtable,  maxdoom);
    removeImpliedSpecNodes(tree, addedNodes);

    return p;
}



// returns the probability of 1 gene giving rise to ngenes after time 'time'
float birthDeathCount(int ngenes, float time, float birthRate, float deathRate)
{
    const float l = birthRate;
    const float u = deathRate;
    const float r = l - u;
    const float a = u / l;

    const float ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time));
    const float p0 = a*ut;
    
    if (ngenes == 0)
        return p0;
    
    // (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
    return (1.0 - p0)*(1.0 - ut) * ipow(ut, ngenes-1);
}




//=============================================================================
// sampling


void getNodeTimes_helper(Node *node, float time, float *times)
{
    const float t = time + node->dist;
    times[node->name] = t;
    
    for (int i=0; i<node->nchildren; i++)
        getNodeTimes_helper(node->children[i], t, times);
}

// gets the time from the root for each node
void getNodeTimes(Tree *tree, float *times)
{
    getNodeTimes_helper(tree->root, 0.0, times);
}

void setNodeTimes(Tree *tree, float *times)
{
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = tree->nodes[i];
        if (node->parent) {
            node->dist = times[node->name] - times[node->parent->name];
        } else {
            // root branch
            node->dist = times[node->name];
        }
    }
}



//  Probability density for for next birth at time 't' given
//  'n' lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
float birthWaitTime(float t, int n, float T, float birth, float death)
{    
    const float l = birth;
    const float u = death;
    const float r = l - u;
    const float a = u / l;

    return n * r * exp(-n*r*t) * \
           pow(1.0 - a * exp(-r * (T - t)), n-1) / \
	   pow(1.0 - a * exp(-r * T), n);
}

//  Probability density for for next birth at time 't' given
//  'n' lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
float birthWaitTime_part(float t, int n, float T, float birth, float death,
		    float denom)
{    
    const float l = birth;
    const float u = death;
    const float r = l - u;
    const float a = u / l;

    return n * r * exp(-n*r*t) * \
           pow(1.0 - a * exp(-r * (T - t)), n-1) / denom;
}

//  Probability density for for next birth at time 't' given
//  'n' lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
float birthWaitTimeDenom(int n, float T, float birth, float death)
{    
    const float l = birth;
    const float u = death;
    const float r = l - u;
    const float a = u / l;

    return pow(1.0 - a * exp(-r * T), n);
}


// Probability of no birth from 'n' lineages starting at time 0, 
// evolving until time 'T' with 'birth' and 'death' rates
// for a reconstructed process.
float probNoBirth(int n, float T, float birth, float death) 
{
    const float l = birth;
    const float u = death;
    const float r = l - u;

    return (1.0 - (l*(1.0 - exp(-r * T)))) / \
	   pow(l - u * exp(-r * T), n);
}


// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n' lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
float sampleBirthWaitTime(int n, float T, float birth, float death)
{
    
    // TODO: could make this much more efficient (use straight line instead of
    // flat line).
    
    // uses rejection sampling
    float denom = birthWaitTimeDenom(n, T, birth, death);
    float start_y = birthWaitTime_part(0, n, T, birth, death, denom);
    float end_y = birthWaitTime_part(T, n, T, birth, death, denom);
    float M = max(start_y, end_y);
    
    while (true) {
        float t = frand(T);
        float f = birthWaitTime_part(t, n, T, birth, death, denom);

        if (frand() <= f / M)
            return t;
    }
}




//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
float birthWaitTime1(float t, float T, float birth, float death,
		    float denom)
{    
    const float l = birth;
    const float u = death;
    const float r = l - u;

    return r * exp(-r*t) / denom;
}

//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
float birthWaitTimeDenom1(float T, float birth, float death)
{    
    const float l = birth;
    const float u = death;
    const float r = l - u;
    const float a = u / l;

    return 1.0 - a * exp(-r * T);
}


// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n'=1 lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
float sampleBirthWaitTime1(float T, float birth, float death)
{
    
    // TODO: could make this much more efficient (use straight line instead of
    // flat line).
    
    // uses rejection sampling
    float denom = birthWaitTimeDenom1(T, birth, death);
    float start_y = birthWaitTime1(0, T, birth, death, denom);
    float end_y = birthWaitTime1(T, T, birth, death, denom);
    float M = max(start_y, end_y);
    
    while (true) {
        float t = frand(T);
        float f = birthWaitTime1(t, T, birth, death, denom);

        if (frand() <= f / M)
            return t;
    }
}






void sampleDupTimes_helper(Node *node, Tree *stree, 
                           int *recon, int *events,
                           float *times, float *stimes,
                           float birthRate, float deathRate)
{     
    
    if (events[node->name] != EVENT_DUP) {
        // speciations happend exactly at species time
        times[node->name] = stimes[recon[node->name]];
    } else {
        if (node->parent == NULL) {
            // TODO: use something realistic
            times[node->name] = stimes[recon[node->name]] - 
                                expovariate(birthRate);
        } else {
            const Node *snode = stree->nodes[recon[node->name]];
            const float start = times[node->parent->name];
            const float difftime = stimes[snode->name] - start;
            
	    // TODO: should not use uniform birth time
            // sample a duplication time
            times[node->name] = start + frand(difftime);
        }
    }

    // recurse
    for (int i=0; i<node->nchildren; i++)
        sampleDupTimes_helper(node->children[i], stree, recon, events, 
			      times, stimes,
			      birthRate, deathRate);
}


void sampleDupTimes(Tree *tree, Tree *stree, int *recon, int *events,
                    float birthRate, float deathRate)
{
    // get species tree times
    float stimes[stree->nnodes];
    getNodeTimes(stree, stimes);
    
    // set gene tree times
    float times[tree->nnodes];
    sampleDupTimes_helper(tree->root, stree, recon, events, times, stimes, 
			  birthRate, deathRate);
    
    // use the times to set branch lengths
    setNodeTimes(tree, times);
}





//=============================================================================
// OLD CODE

// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
float birthDeathCount2(int start, int end, float time, float birthRate, float deathRate)
{
    const float l = birthRate;
    const float u = deathRate;
    const float r = l - u;
    const float a = u / l;

    const float ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time));
    const float p0 = a*ut;
    
    // all 'start' genes die out
    if (end == 0) {
        return ipow(p0, start);
    }
    
    const int iter = (start < end) ? start : end;
    float p = 0.0;
    for (int j=0; j<iter; j++) {
        p += choose(start, j) *
             choose(start + end - j - 1, start - 1) *
             ipow(p0, start-j) *
             ipow(ut, end-j) *
             ipow(1 - p0 - ut, j);
    }
        
    return p;
}




float birthDeathDensityNoExtinct(float *times, int ntimes, float maxtime, 
                                 float birthRate, float deathRate)
{
    const float l = birthRate;
    const float u = deathRate;
    const float r = l - u;
    const float a = u / l;
    const float T = maxtime;
            
    assert(ntimes >= 2);

    float p = 1.0;

    // (n - 1)! r^(n-2) (1 - a)^n
    for (int i=1; i<=ntimes-2; i++)
        p *= i * r * (1.0 - a);
    p *= (ntimes - 1) * (1.0 - a) * (1.0 - a);

    // exp( r * sum_{i=3}^n x_i)
    float sum = 0.0;
    for (int i=2; i<ntimes; i++)
        sum += T - times[i];
    p *= exp(r * sum);

    // prod_{i=2}^n 1/(exp(r x_i) - a)^2
    float prod = 1.0;
    for (int i=1; i<ntimes; i++) {
        float tmp = exp(r*(T-times[i])) - a;
        prod *= (tmp*tmp);
    }
    p *= 1.0 / prod;
    
    return p;
}


float birthDeathDensity(float *times, int ntimes, float maxtime, 
                        float birthRate, float deathRate)
{
    const float l = birthRate;
    const float u = deathRate;
    const float r = l - u;
    const float a = u / l;
    const float T = maxtime;

    const float uT = (1.0 - exp(-r*T)) / (1.0 - a * exp(-r*T));
    const float p0 = a*uT;
    
    
    
    if (ntimes == 0) {
        return p0;
    } else if (ntimes == 1) {
        const float p1 = (1.0 - a*uT) * (1.0 - uT);       
        return p1;
    } else {
        float p = (1.0 - p0);
        
        // (n - 1)! r^(n-2) (1 - a)^n
        for (int i=1; i<=ntimes-2; i++)
            p *= i * r * (1.0 - a);
        p *= (ntimes - 1) * (1.0 - a) * (1.0 - a);
        
        //printf("p_1 %f\n", p);
        
        // exp( r * sum_{i=3}^n x_i)
        float sum = 0.0;
        for (int i=2; i<ntimes; i++)
            sum += T - times[i];
        p *= exp(r * sum);
        
        // prod_{i=2}^n 1/(exp(r x_i) - a)^2
        float prod = 1.0;
        for (int i=1; i<ntimes; i++) {
            float tmp = exp(r*(T-times[i])) - a;
            prod *= (tmp*tmp);
        }
        p *= 1.0 / prod;

        // f(t_2 - t_1, 1)
        p *= (r * exp(-r * (times[1] - times[0]))) / 
             (1.0 - a * exp(-r * (T - times[0])));
        
        for (int i=1; i<=ntimes-1; i++)        
            p *= 2.0 / (i * i);
        p /= ntimes;
        
        
        //printFloatArray(times, ntimes);
        //printf("p %f; T %f\n", p, T);
        return p;
    }
}


/*

void birthDeathTreePrior_recurse(Node *node, float time, int *events, 
                                 ExtendArray<float> &times)
{
    if (events[node->name] == EVENT_DUP) {
        times.append(time += node->dist);
        
        for (int j=0; j<node->nchildren; j++) {
            birthDeathTreePrior_recurse(node->children[j], time + node->dist,
                                        events, times);
        }
    }
}




// NOTE: assumes binary species tree
float birthDeathTreePrior(Tree *tree, SpeciesTree *stree, int *recon, 
                          int *events, float birthRate, float deathRate)
{
    float prob = 0.0;
    float _times[tree->nnodes];
    ExtendArray<float> times(0, tree->nnodes, _times);
    
    // catch undefined params
    if (birthRate == deathRate)
        deathRate = .99 * birthRate;    
    
    // preroot duplications
    if (events[tree->root->name] == EVENT_DUP) {
        times.clear();
        times.append(0.0);
        birthDeathTreePrior_recurse(tree->root, 0.0, events, times);
        float maxtime = 0.0;
        for (int i=0; i<times.size(); i++)
            if (times[i] > maxtime)
                maxtime = times[i];
        prob += log(birthDeathDensityNoExtinct(times, times.size(),
                                               maxtime, birthRate, deathRate));
        //int rootdups = countRootDups(tree->root, events);
        //prob += log(birthDeathCount(rootdups+1, pretime, birthRate, deathRate));
    }
 

    // loop through speciation nodes in tree2
    for (int i=0; i<tree->nnodes; i++) {
        const Node *node = tree->nodes[i];
        
        if (events[node->name] == EVENT_SPEC) {
            if (node->nchildren == 1) {
                // calc loss prob
                const Node *snode = stree->nodes[recon[node->name]];
                float maxtime;
                
                if (snode->children[0]->name == 
                    recon[node->children[0]->name])
                {
                    maxtime = snode->children[1]->dist;
                } else {
                    maxtime = snode->children[0]->dist;
                }
                prob += log(birthDeathDensity(NULL, 0, maxtime, 
                                              birthRate, deathRate));
                
            }
            
            for (int j=0; j<node->nchildren; j++) {
                const float maxtime =
                    stree->nodes[recon[node->children[j]->name]]->dist;
            
                if (events[node->children[j]->name] != EVENT_DUP) {
                    // no apparent dup/loss  1->1
                    prob += log(birthDeathDensity(NULL, 1, maxtime, 
                                                  birthRate, deathRate));
                } else {
                    // duplication
                    times.clear();
                    times.append(0.0);
                    
                    birthDeathTreePrior_recurse(node->children[j], 0.0, events, 
                                                times);
                    
                    //printf("maxtime %f\n", maxtime);
                    
                    prob += log(birthDeathDensity(times, times.size(),
                                                  maxtime, birthRate, deathRate));
                }
            }
        }
    }
    
    
    // clean up
    times.detach();
    
    return prob;
}


//=============================================================================
// DEAD CODE: incorrectly hashed redundantly labeled trees 

class HashNode2 {
public:
    HashNode2(int nodeid, int start, int end, int *key) :
	nodeid(nodeid),
	start(start),
	end(end),
	len(end - start),
	key(key)
    {
    }
    int nodeid;
    int start;
    int end;
    int len;
    int *key;
};


int hashNodeCmp2(const void *_a, const void *_b)
{
    HashNode2 *a = *((HashNode2**) _a);
    HashNode2 *b = *((HashNode2**) _b);
    
    // first compare diffs
    int diff = a->len - b->len;
    if (diff) {
	return diff;
    } else {
	int *keya = a->key + a->start;
	int *keyb = b->key + b->start;
	for (int i=0; i <= a->len; i++) {
	    diff = keya[i] - keyb[i];
	    if (diff)
		return diff;
	}

	// both hashes are the same
	// 1. same length
	// 2. same key subsequence
	return 0;
    }
}



void getHashIds2(Tree *tree, int *recon, int *hashids)
{
    // get post order of nodes
    ExtendArray<Node*> postnodes(0, tree->nnodes);
    getTreePostOrder(tree, &postnodes);
    
    // order children
    ExtendArray<int> ordering(tree->nnodes);
    for (int i=0; i<postnodes.size(); i++)
    {
        Node *node=postnodes[i];
        
        if (node->isLeaf()) {
            ordering[node->name] = recon[node->name];
        } else {
            // propogate the min order to the parent
            int minorder = ordering[node->children[0]->name];
            for (int j=1; j<node->nchildren; j++) {
                int order = ordering[node->children[j]->name];
                if (order < minorder)
                    minorder = order;
            }
            ordering[node->name] = minorder;
        }
    }
    
    // get a sorted post ordering of nodes
    ExtendArray<Node*> sortpostnodes(0, tree->nnodes);
    getTreeSortedPostOrder(tree, &sortpostnodes, ordering);
    
    // generate a unique key for this topology
    // postfix notation for a tree
    // ((A,B),(C)) is represented as
    // A, B, -1, C, -2, -1
    ExtendArray<int> key(tree->nnodes);
    ExtendArray<int> start_stack(0, tree->nnodes);
    ExtendArray<HashNode2*> hashnodes(0, tree->nnodes);

    for (int i=0; i<sortpostnodes.size(); i++) {
        Node *node = sortpostnodes[i];
	int last;
        if (node->isLeaf()) {
            key[i] = recon[node->name];
	    last = i;
        } else if (node->nchildren == 2) {
            key[i] = -1;
	    start_stack.pop();
	    last = start_stack.pop();
	} else if (node->nchildren == 1) {
	    key[i] = -2;
	    last = start_stack.pop();
	} else {
	    // cannot handle trifurcating trees
	    assert(0);
	}

	hashnodes.append(new HashNode2(node->name, last, i, key));
	start_stack.append(last);
    }

    //printf("key = ");
    //printIntArray(key, key.size());
    //printf("\n");
    
    // sort hashnodes
    qsort((void*) hashnodes.get(), hashnodes.size(), sizeof(HashNode2*),
	  hashNodeCmp2);    
    
    int hashid = 0;
    hashids[hashnodes[0]->nodeid] = hashid;
    for (int i=1; i<hashnodes.size(); i++) {
	// use new hashid if nodes differ
	if (hashNodeCmp2(&hashnodes[i], &hashnodes[i-1]))
	    hashid++;
	hashids[hashnodes[i]->nodeid] = hashid;
    }
    

    // clean up
    for (int i=0; i<hashnodes.size(); i++)
	delete hashnodes[i];

    //printf("hashids = ");
    //printIntArray(hashids, tree->nnodes);
    //printf("\n");
}

//=============================================================================

*/


// NOTE: assumes binary species tree
// branch lengths are ignored, only counts are used
float birthDeathTreeQuickPrior(Tree *tree, SpeciesTree *stree, int *recon, 
                               int *events, float birthRate, float deathRate)
{   

    float prob = 0.0;
    
    // catch undefined params
    if (birthRate == deathRate)
        deathRate = .99 * birthRate;    
    
    // preroot duplications
    if (events[tree->root->name] == EVENT_DUP) {
        // ignore root right now
        const float predupprob = 0.02;
        int rootdups = countDups(tree->root, events);
        
        prob += log(predupprob) * rootdups + log(1.0 - predupprob);
        //prob += log(birthDeathCount(rootdups+1, pretime, birthRate, deathRate));
    }
 

    // loop through speciation nodes in tree2
    for (int i=0; i<tree->nnodes; i++) {
        const Node *node = tree->nodes[i];
        
        if (events[node->name] == EVENT_SPEC) {
            if (node->nchildren == 1) {
                // calc loss prob
                const Node *snode = stree->nodes[recon[node->name]];
                float maxtime;
                
                if (snode->children[0]->name == 
                    recon[node->children[0]->name])
                    maxtime = snode->children[1]->dist; 
                else
                    maxtime = snode->children[0]->dist;
                assert(maxtime > 0.0);
                prob += log(birthDeathCount(0, maxtime, birthRate, deathRate));
            }
            
            for (int j=0; j<node->nchildren; j++) {
                const float maxtime =
                    stree->nodes[recon[node->children[j]->name]]->dist;
                assert(maxtime > 0.0);
            
                if (events[node->children[j]->name] != EVENT_DUP) {
                    // no apparent dup/loss  1->1
                    prob += log(birthDeathCount(1, maxtime, 
                                                birthRate, deathRate));
                } else {
                    // duplication
                    int count = countDups(node->children[j], events);
                    
                    //printf("maxtime %f\n", maxtime);
                    
                    prob += log(birthDeathCount(count+1, maxtime, 
                                                birthRate, deathRate));
                }
            }
        }
    }
    
    return prob;
}



}

} // namespace spidir
