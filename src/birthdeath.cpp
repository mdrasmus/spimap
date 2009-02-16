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


// TODO: does not handle branches above the species tree root yet
// NOTE: assumes binary species tree
float birthDeathTreePrior(Tree *tree, Tree *stree, int *recon, 
                          int *events, float birthRate, float deathRate,
                          float *doomtable, int maxdoom)
{

    float prob = 0.0;
    ExtendArray<Node*> subleaves(0, tree->nnodes);
    
    // catch undefined params
    if (birthRate == deathRate)
        deathRate = .99 * birthRate;
    
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

                float nhist = numSubtopologyHistories(tree, node, subleaves);
                float thist = numHistories(subleaves.size());

                // sum over ndoom
                float sum = 0.0;
                for (int ndoom=0;  ndoom<=maxdoom; ndoom++) {
                    sum += (birthDeathCount(subleaves.size() + ndoom, 
                                            schild->dist,
                                            birthRate, deathRate) *
                            ipow(exp(doomtable[schild->name]), ndoom));
                }

                prob += log(nhist) - log(thist) + log(sum);
            }
        }
    }
    
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
    printf("> %d\n", addedNodes);
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
