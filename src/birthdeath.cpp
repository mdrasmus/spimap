
// c/c++ includes
#include <math.h>

// spidir includes
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


/*
int numRedundantTopologies_helper(Node *node, 
				 bool *leafset, 
				 int *hashids)
{
    if (leafset[node->name] == true)
	return 0;
    else {
	if (node->nchildren == 1)
	    return numRedundantTopologies_helper(node->children[0],
						leafset, hashids);
	else if (node->nchildren == 2) {
	    int n = int(hashids[node->children[0]->name] == 
			hashids[node->children[1]->name]);
	    return n + numRedundantTopologies_helper(node->children[0],
						    leafset, hashids)
		     + numRedundantTopologies_helper(node->children[1],
						    leafset, hashids);
	} else {
	    // NOTE: cannot handle multifurcating nodes
	    assert(0);
	}
    }
}
*/

// N_n(T, allLeaves=false) = 2^{-M(T)} prod_u c_u!
// N_n(T, allLeaves=true) = 2^{-M(T)} L(T)!
double numRedundantTopologies(Tree *tree, Node *root, 
                              ExtendArray<Node*> &leaves, 
                              int *hashids, bool allLeaves)
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

    double val = 0.0;

    if (allLeaves) {
        // val = log(factorial(leaves.size()))
        for (double i=2.0; i<=leaves.size(); i+=1.0)
            val += logf(i);
    } else {
        // get hashes
        ExtendArray<int> leafhashes(0, leaves.size());
        for (int i=0; i<leaves.size(); i++)
            leafhashes.append(hashids[leaves[i]->name]);
    
        qsort((void*) leafhashes.get(), leafhashes.size(), sizeof(int), intCmp);

        double colorsize = 1;
        for (int i=1; i<leaves.size(); i++) {
            if (leafhashes[i] != leafhashes[i-1]) {
                // val *= factorial(colorsize)
                for (double j=2; j<=colorsize; j+=1.0)
                    val += logf(j);
                colorsize = 1.0;
            } else {
                colorsize += 1.0;
            }
        }
        for (double j=2; j<=colorsize; j+=1.0)
            val += logf(j);
    }

    // divide by 2^M
    for (int i=0; i<nmirrors; i++)
        val -=  logf(2.0);

    return val;
}


/*
// N_c = multinomial{|L(T)|!}{c_1 ... c_k}
double numRedundantTopologies2(Tree *tree, Node *root, 
			     ExtendArray<Node*> &leaves, 
			     int *hashids)
{

    // get hashes
    ExtendArray<int> leafhashes(0, leaves.size());
    for (int i=0; i<leaves.size(); i++)
	leafhashes.append(hashids[leaves[i]->name]);
    
    qsort((void*) leafhashes.get(), leafhashes.size(), sizeof(int), intCmp);

    // prod_u c_u!
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

    // L(T)!
    double num = 1.0;
    for (double j=2; j<=leaves.size(); j+=1.0)
	num *= j;

    return num / val;
}
*/


double birthDeathTopology(Node *node, float birthRate, float deathRate,
                         ExtendArray<Node*> &leaves)
{
    // numTopologyHistories(subtree) / numHistories(ngenes)

    return numHistories(leaves.size());
}

// computes the entries of the doom probabilty table
void calcDoomTable(Tree *tree, float birthRate, float deathRate, int maxdoom,
                   double *doomtable)
{
    // get nodes in post order
    ExtendArray<Node*> nodes(0, tree->nnodes);
    getTreePostOrder(tree, &nodes);

    
    for (int i=0; i<tree->nnodes; i++) {
        Node *node = nodes[i];
        
        if (node->isLeaf()) {
            doomtable[node->name] = -INFINITY;
        } else {
            double prod = 0.0;
            
            for (int j=0; j<node->nchildren; j++) {
                Node *child = node->children[j];
                double sum = 0.0;

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
double birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                           int *events, float birthRate, float deathRate,
                           double *doomtable, int maxdoom)
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
		    nhist = numSubtopologyHistories(tree, node, subleaves);
		    thist = numHistories(subleaves.size());		    

		    if (subleaves[0]->isLeaf()) {
                        // correct subtrees that have leaves
			nhist *= exp(numRedundantTopologies(tree, node, 
                                                            subleaves, 
                                                            hashids,
                                                            false));
		    } else {
                        // correct subtrees that have leaves
			double a = exp(numRedundantTopologies(tree, node, 
                                                              subleaves, 
                                                              hashids,
                                                              true));
                        nhist *= a;
                    }
		}
		
                // sum over ndoom
                double sum = 0.0;
                for (int ndoom=0;  ndoom<=maxdoom; ndoom++) {
                    int totleaves = subleaves.size() + ndoom;
                    sum += fchoose(totleaves, ndoom) *
                        birthDeathCount(totleaves, 
                                        schild->dist,
                                        birthRate, deathRate) *
                        ipow(exp(doomtable[schild->name]), ndoom);
                }

                prob += log(nhist) - log(thist) + log(sum);
            }
        }
    }

    ExtendArray<Node*> leaves(0, tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
	if (tree->nodes[i]->isLeaf())
	    leaves.append(tree->nodes[i]);
    }
    double x = numRedundantTopologies(tree, tree->root, leaves, 
                                      hashids, false);
    prob -= x;
    
    return prob;
}



// Convenience function
// Adds and removes implied species nodes to the gene tree
// NOTE: assumes binary species tree
double birthDeathTreePriorFull(Tree *tree, Tree *stree, int *recon, 
                              int *events, float birthRate, float deathRate,
                              double *doomtable, int maxdoom)
{
    ExtendArray<int> recon2(0, 2 * tree->nnodes);
    recon2.extend(recon, tree->nnodes);

    ExtendArray<int> events2(0, 2 * tree->nnodes);
    events2.extend(events, tree->nnodes);


    int addedNodes = addImpliedSpecNodes(tree, stree, recon2, events2);
    double p = birthDeathTreePrior(tree, stree, recon2, events2, 
                                   birthRate, deathRate,
                                   doomtable,  maxdoom);
    removeImpliedSpecNodes(tree, addedNodes);

    return p;
}



// returns the probability of 1 gene giving rise to ngenes after time 'time'
double birthDeathCount(int ngenes, float time, 
                       float birthRate, float deathRate)
{
    const double l = birthRate;
    const double u = deathRate;

    if (birthRate == deathRate) {
        const double ut = time / (1.0 / birthRate + time);
        if (ngenes == 0)
            return ut;
        return ((1.0 - ut)*(1.0 - ut)) * ipow(ut, ngenes-1);
    }

    const double r = l - u;
    const double a = u / l;

    const double ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time));
    const double p0 = a*ut;

    if (ngenes == 0)
        return p0;
    
    // (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
    return (1.0 - p0)*(1.0 - ut) * ipow(ut, ngenes-1);
}


// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
// slower more stable computation
double birthDeathCounts2(int start, int end, float time, 
                         float birth, float death)
{
    const double l = birth;
    const double u = death;
    const double r = l - u;
    const double a = u / l;

    const double ertime = exp(-r*time);
    const double ut = (1.0 - ertime) / (1.0 - a * ertime);
    const double p0 = a*ut;
    
    // all 'start' genes die out
    if (end == 0) {
        return ipow(p0, start);
    }
    
    const int iter = (start < end) ? start : end;
    double p = 0.0;

    for (int j=0; j<=iter; j++) {
        //printf("j %d\n", j);
        p += fchoose(start, j) *
             fchoose(start + end - j - 1, start - 1) *
             ipow(p0, start-j) *
             ipow(ut, end-j) *
             ipow(1 - p0 - ut, j);
    }
    
    // do not allow invalid values to propogate
    if (isnan(p) || isinf(p) || p > 1.0) {
        printf("p=%e genes=(%d, %d) b=%f d=%f t=%f\n", p, start, end,
               birth, death, time);
        fflush(stdout);
        assert(0);
    }
    return p;
}

// returns the probability of 'start' genes giving rise to 'end' genes after 
// time 'time'
// much faster computation than birthDeathCounts2
double birthDeathCounts(int start, int end, float time, 
                       float birth, float death)
{
    if (start == 0) {
        if (end == 0)
            return 1.0;
        else
            return 0.0;
    }

    const double ertime = exp((birth-death)*time);
    const double tmp = (ertime-1.0) / (birth*ertime - death);
    const double a = death * tmp;
    const double b = birth * tmp;
    
    // all 'start' genes die out
    if (end == 0) {
        return ipow(a, start);
    }
    
    // compute base case
    double f = ipow(a, start) * ipow(b, end);
    if (start > 1)
        f *= (start + end - 1);
    for (int k=2; k<start; k++)
        f *= (start + end - k) / double(k);


    double p = f;
    double x = start;
    double y = end;
    double z = start + end - 1;
    const double oneab = 1.0 - a - b;
    const int iter = (start < end) ? start : end;
    for (int j=1; j<=iter; j++) {
        f *= (oneab * x * y / (j * a * b * z));
        //printf("p=%e f=%e j=%d\n", p, f, j);
        p += f;
        x--;
        y--;
        z--;
    }

    if (p < 0.0)
        p = 0.0;
    

    if (p > 1.0) {
        // resort to a slower more stable function
        return birthDeathCounts2(start, end, time, birth, death);
    }

    return p;
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
double birthWaitTime(float t, int n, float T, float birth, float death)
{    
    if (birth == death) {
        const double t2 = t - T;
        const double nl = 1.0 / birth;
        return birth * n * ipow(-nl+t2, n) / ipow(-nl-T, n) / (1.0-birth*t2);
    }

    const double r = birth - death;
    const double a = death / birth;

    return n * r * exp(-n*r*t) * \
           pow(1.0 - a * exp(-r * (T - t)), n-1) / \
	   pow(1.0 - a * exp(-r * T), n);
}

// numerator for birthWaitTime
double birthWaitTimeNumerator(float t, int n, float T, 
                              float birth, float death,
                              double denom)
{    
    const double r = birth - death;
    const double a = death / birth;

    return n * r * exp(-n*r*t) * 
           pow(1.0 - a * exp(-r * (T - t)), n-1) / denom;
}

// denominator for birthWaitTime
double birthWaitTimeDenom(int n, float T, float birth, float death)
{    
    const double r = birth - death;
    const double a = death / birth;

    return pow(1.0 - a * exp(-r * T), n);
}


// Probability of no birth from 'n' lineages starting at time 0, 
// evolving until time 'T' with 'birth' and 'death' rates
// for a reconstructed process.
double probNoBirth(int n, float T, float birth, float death) 
{
    if (birth == 0.0)
        return 1.0;
    else if (birth == death) {
        return 1.0 / ipow(1.0 + birth * T, n);
    }

    const double r = birth - death;
    const double exprt = exp(-r * T);

    return pow(1.0 - (birth*(1.0 - exprt)) / (birth - death * exprt), n);
}



// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n' lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
double sampleBirthWaitTime(int n, float T, float birth, float death)
{
    
    // TODO: could make this more efficient
    
    if (birth == death) {
        double start_y = birthWaitTime(0, n, T, birth, death);
        double end_y = birthWaitTime(T, n, T, birth, death);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTime(t, n, T, birth, death);

            if (frand() <= f / M)
                return t;
        }
    } else {
        // uses rejection sampling
        double denom = birthWaitTimeDenom(n, T, birth, death);
        double start_y = birthWaitTimeNumerator(0, n, T, birth, death, denom);
        double end_y = birthWaitTimeNumerator(T, n, T, birth, death, denom);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTimeNumerator(t, n, T, birth, death, denom);

            if (frand() <= f / M)
                return t;
        }
    }
 }




//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTime1(float t, float T, float birth, float death,
		    float denom)
{
    const double r = birth - death;

    return r * exp(-r*t) / denom;
}

//  Probability density for for next birth at time 't' given
//  'n'=1 lineages starting at time 0, evolving until time 'T' with a
//  'birth' and 'death' rates for a reconstructed process.
double birthWaitTimeDenom1(float T, float birth, float death)
{
    const double r = birth - death;
    const double a = death / birth;

    return 1.0 - a * exp(-r * T);
}


// Sample the next birth event from a reconstructed birthdeath process.
// Let there be 'n'=1 lineages at time 0 that evolve until time 'T' with
// 'birth' and 'death' rates.
// Conditioned that a birth will occur
double sampleBirthWaitTime1(float T, float birth, float death)
{    
    // TODO: could make this much more efficient

    if (birth == death) {
        // uses rejection sampling
        double start_y = birthWaitTime(0, 1, T, birth, death);
        double end_y = birthWaitTime(T, 1, T, birth, death);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTime(t, 1, T, birth, death);

            if (frand() <= f / M)
                return t;
        }

    } else {
        // uses rejection sampling
        double denom = birthWaitTimeDenom1(T, birth, death);
        double start_y = birthWaitTime1(0, T, birth, death, denom);
        double end_y = birthWaitTime1(T, T, birth, death, denom);
        double M = max(start_y, end_y);
    
        while (true) {
            double t = frand(T);
            double f = birthWaitTime1(t, T, birth, death, denom);

            if (frand() <= f / M)
                return t;
        }
    }
}







}

} // namespace spidir
