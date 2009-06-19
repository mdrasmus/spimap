//=============================================================================
//  SPIDIR - tree search


// 3rd party includes
#include <algorithm>


// spidir includes
#include "common.h"
#include "branch_prior.h"
#include "parsimony.h"
#include "phylogeny.h"
#include "search.h"
#include "seq_likelihood.h"
#include "Matrix.h"


namespace spidir {

//=============================================================================
// Nearest Neighbor Interchange Topology Proposal

/*

    Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      nodeb    node1
               /  \
         nodea     * 

*/
void performNni(Tree *tree, Node *nodea, Node *nodeb)
{
    Node *node1 = nodea->parent;
    Node *node2 = nodeb->parent;
    
    // assert that node1 and node2 are incident to the same branch
    assert(node1->parent == node2 ||
           node2->parent == node1);
    
    // find child indexes
    int a = (node1->children[0] == nodea) ? 0 : 1;
    assert(node1->children[a] == nodea);

    int b = (node2->children[0] == nodeb) ? 0 : 1;
    assert(node2->children[b] == nodeb);
    
    // swap parent pointers
    nodea->parent = node2;
    nodeb->parent = node1;
    
    // swap child pointers
    node2->children[b] = nodea;
    node1->children[a] = nodeb;
}


void proposeRandomNni(Tree *tree, Node **a, Node **b)
{
    // find edges for NNI
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    Node *node1 = tree->nodes[choice];
    Node *node2 = tree->nodes[choice]->parent;
    *a = node1->children[irand(2)];
    *b = (node2->children[0] == node1) ? node2->children[1] :
                                         node2->children[0];
}


//=============================================================================
// Subtree pruning and regrafting (SPR)


/*
    a = subtree
    e = newpos

    BEFORE
            ....
        f         d
       /           \
      c             e
     / \           ...
    a   b
   ... ...

    AFTER

        f         d
       /           \
      b             c
     ...           / \
                  a   e
                 ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary

*/
void performSpr(Tree *tree, Node *subtree, Node *newpos)
{
    Node *a = subtree;
    Node *e = newpos;

    Node *c = a->parent;
    Node *f = c->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    const int ci = (f->children[0] == c) ? 0 : 1;
    Node *d = e->parent;
    const int ei = (d->children[0] == e) ? 0 : 1;

    d->children[ei] = c;
    c->children[bi] = e;
    f->children[ci] = b;
    b->parent = f;
    b->dist += c->dist;
    e->dist /= 2.0;
    c->dist = e->dist;
    c->parent = d;
    e->parent = c;
}

/*
    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
    
          d
         / \
        e  ...
       / \
      c  ...         
     / \           
    a   b
   ... ...

    AFTER
          d
         / \
        c
       / \
      a   e
     ... / \
        b  ...
       ...
       
  What if d == f  (also equivalent to NNI) this is OK
  
    BEFORE
          
        f
       / \
      c   e
     / \  ...
    a   b
   ... ...

    AFTER
          
        f
       / \
      b   c  
     ... / \ 
        a   e
       ... ...  
*/



/*
    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary
*/
void proposeRandomSpr(Tree *tree, Node **subtree, Node **newpos)
{
    assert(tree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    *subtree = a;
    
    // find sibling (b) of a
    Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    
    // choose newpos (e)
    Node *e = NULL;
    do {
        choice = irand(tree->nnodes);
        e = tree->nodes[choice];
        
        // test if e is a valid choice
        if (e->parent == NULL || e == a || e == c || e == b)
            continue;
        
        // also test if e is a descendent of a
        bool under_a = false;
        for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent) {
            if (ptr == a) {
                under_a = true;
                break;
            }
        }            
        
        if (under_a)
            continue;
        
        break;
    } while (true);
    *newpos = e;
}



//=============================================================================
// NNI Proposer

NniProposer::NniProposer(SpeciesTree *stree, int *gene2species,
                         int niter) :
    nodea(NULL),
    nodeb(NULL),
    nodec(NULL),
    noded(NULL),
    oldroot1(NULL),
    oldroot2(NULL),
    stree(stree),
    gene2species(gene2species),
    niter(niter),
    iter(0)
{}
    
    
void NniProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;
    const float doubleNniProb = 0.0; // NO double NNI
    
    // increase iteration
    iter++;
    
    nodea = nodeb = nodec = noded = NULL;

    // propose new tree
    proposeRandomNni(tree, &nodea, &nodeb);
    performNni(tree, nodea, nodeb);

    if (frand() < doubleNniProb) {
        //printLog(LOG_MEDIUM, "search: double NNI\n");
        proposeRandomNni(tree, &nodec, &noded);
        performNni(tree, nodec, noded);
    }
    
    // reroot tree if stree is given
    if (frand() < rerootProb) {
        oldroot1 = tree->root->children[0];
        oldroot2 = tree->root->children[1];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        }
    } else {
        oldroot1 = NULL;
        oldroot2 = NULL;
    }
    
    // debug: keep track of correct tree in search
    testCorrect(tree);
    
    //assert(tree->assertTree());
}

void NniProposer::revert(Tree *tree)
{
    // reject, undo topology change


    if (oldroot1)
        tree->reroot(oldroot1, oldroot2);
    
    if (nodec)
        performNni(tree, nodec, noded);
    if (nodea)
        performNni(tree, nodea, nodeb);
}


bool NniProposer::more()
{
    return iter < niter;
}



//=============================================================================
// SPR Proposer

SprProposer::SprProposer(SpeciesTree *stree, int *gene2species,
                         int niter) :
    NniProposer(stree, gene2species, niter)
{
}
    
    
void SprProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;

    // increase iteration
    iter++;
    
    // choose a SPR move
    proposeRandomSpr(tree, &nodea, &nodeb);
    
    // remember sibling of nodea
    const Node *p = nodea->parent;
    nodec = (p->children[0] == nodea) ? p->children[1] : p->children[0];
    
    // perform SPR move
    performSpr(tree, nodea, nodeb);

    // reroot tree if stree is given
    if (frand() < rerootProb) {
        oldroot1 = tree->root->children[0];
        oldroot2 = tree->root->children[1];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        }
    } else {
        oldroot1 = NULL;
        oldroot2 = NULL;
    }

    // debug: keep track of correct tree in search
    testCorrect(tree);
}

void SprProposer::revert(Tree *tree)
{

    if (oldroot1)
        tree->reroot(oldroot1, oldroot2);

    performSpr(tree, nodea, nodec);
}


//=============================================================================
// SPR + NNI Proposer

SprNniProposer::SprNniProposer(SpeciesTree *stree, int *gene2species,
                               int niter, float sprRatio) :
    NniProposer(stree, gene2species, niter),
    sprRatio(sprRatio)
{
}
    
    
void SprNniProposer::propose(Tree *tree)
{
    const float noChange = 0.00;
    float choice = frand();
    
    if (correctTree) {
        if (tree->sameTopology(correctTree))
            correctSeen = true;
    }
    
    if (choice < noChange) {
        lastPropose = PROPOSE_NONE;
        iter++; // increase iteration
        
    } if (choice < noChange + (1.0 - noChange) * sprRatio) {
        lastPropose = PROPOSE_SPR;
        
        // increase iteration
        iter++;
        
        proposeRandomSpr(tree, &nodea, &nodeb);
        
        // remember sibling of nodea
        const Node *p = nodea->parent;
        nodec = (p->children[0] == nodea) ? p->children[1] : p->children[0];
        
        performSpr(tree, nodea, nodeb);
        
    } else {
        lastPropose = PROPOSE_NNI;
        NniProposer::propose(tree);
    }
}

void SprNniProposer::revert(Tree *tree)
{
    if (lastPropose == PROPOSE_NONE)
        return;
    else if (lastPropose == PROPOSE_SPR) {
        performSpr(tree, nodea, nodec);
    
    } else if (lastPropose == PROPOSE_NNI) {
        NniProposer::revert(tree);
    } else {
        assert(false);
    }
}



//=============================================================================
// SPR Neighborhood Proposer

SprNbrProposer::SprNbrProposer(SpeciesTree *stree, int *gene2species,
                               int niter, int radius) :
    NniProposer(stree, gene2species, niter),
    radius(radius),
    basetree(NULL)
{
}
    
    
void SprNbrProposer::propose(Tree *tree)
{
    // ensure the same tree is used for each proposal
    if (!basetree)
        basetree = tree;
    else
        assert(basetree == tree);


    // start a new subtree
    if (iter == 0 || queue.size() == 0) {
        iter = 0;
        pickNewSubtree();
    }

    iter++; // increase iteration
    
    // get new branch point
    nodea = queue.front();
    queue.pop_front();
    
    // remember sibling of subtree (nodeb)
    const Node *p = subtree->parent;
    nodeb = (p->children[0] == subtree) ? p->children[1] : p->children[0];
    
    // perform SPR move
    performSpr(tree, subtree, nodea);
    testCorrect(tree); // debug: keep track of correct tree in search

    tree->assertTree();
}

void SprNbrProposer::revert(Tree *tree)
{
    performSpr(tree, subtree, nodeb);
    tree->assertTree();
}


void SprNbrProposer::pickNewSubtree()
{
    const Tree *tree = basetree;

    assert(basetree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    Node *a = tree->nodes[choice];
    subtree = a;
    
    // find sibling (b) of a
    Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    
    // uninitialize path distances
    pathdists.clear();
    for (int i=0; i<tree->nnodes; i++)
        pathdists.push_back(-1);
    
    // setup path distances and queue
    queue.clear();
    pathdists[a->name] = 0;
    pathdists[c->name] = 0;
    pathdists[b->name] = 0;
    list<Node*> tmpqueue;
    tmpqueue.push_back(c);
    tmpqueue.push_back(b);

    // traverse tree via depth first traversal
    while (tmpqueue.size() > 0) {
        Node *n = tmpqueue.front();
        tmpqueue.pop_front();
        
        // do not traverse beyond radius
        if (pathdists[n->name] >= radius)
            continue;

        // queue only valid new branch points:
        // n must not be root, a, descendant of a, c (parent of a), or  
        // b (sibling of a)
        if (n->parent && n != subtree && n != b && n != c) {
            queue.push_back(n);
        }

        // queue up unvisited neighboring edges
        Node *w = n->parent;

        if (w && pathdists[w->name] == -1) {
            pathdists[w->name] = pathdists[n->name] + 1;
            tmpqueue.push_back(w);
        }

        if (n->nchildren == 2) {
            Node *u = n->children[0];
            Node *v = n->children[1];

            if (pathdists[u->name] == -1) {
                pathdists[u->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(u);
            }

            if (pathdists[v->name] == -1) {
                pathdists[v->name] = pathdists[n->name] + 1;
                tmpqueue.push_back(v);
            }
        }
    }
}


void SprNbrProposer::reset() { iter = 0; }



//=============================================================================
// Dup/Loss proposer

DupLossProposer::DupLossProposer(TopologyProposer *proposer, 
                                 SpeciesTree *stree, int *gene2species,
                                 float dupprob, float lossprob,
                                 int quickiter, int niter,
                                 int nsamples, bool extend) :
    proposer(proposer),
    quickiter(quickiter),
    niter(niter),
    iter(0),
    correctTree(NULL),
    correctSeen(false),
    stree(stree),
    gene2species(gene2species),
    dupprob(dupprob),
    lossprob(lossprob),
    recon(0),
    events(0),
    oldtop(NULL),
    nsamples(nsamples),
    samplei(0),
    treesize(0)
{
    doomtable = new float [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, maxdoom, doomtable);

}


DupLossProposer::~DupLossProposer()
{
    //if (oldtop)
    //    delete oldtop;
    delete [] doomtable;

    //clearQueue();
}

void DupLossProposer::propose(Tree *tree)
{
    iter++;
    
    // record correct tree if we visit it
    if (correctTree) {
        if (tree->sameTopology(correctTree))
            correctSeen = true;
    }
    
    
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->propose(tree);
        return;
    }
    
    // save old topology
    oldtop = tree->copy();
    
    // do simple proposal
    if (frand() < .2) {
        proposer->propose(tree);
        return;
    }
    
    // recon tree to species tree
    recon.ensureSize(tree->nnodes);
    events.ensureSize(tree->nnodes);
    recon.setSize(tree->nnodes);
    events.setSize(tree->nnodes);
    
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    
    ExtendArray<Tree*> trees(0, quickiter);
    ExtendArray<float> logls(0, quickiter);
    
    float sum = -INFINITY;
    
    // make many subproposals
    for (int i=0; i<quickiter; i++) {
        proposer->propose(tree);
        
        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
    
        int addedNodes = addImpliedSpecNodes(tree, stree, recon, events);
        float logl = birthDeathTreeQuickPrior(tree, stree, recon, events, 
                                              dupprob, lossprob);
        sum = logadd(sum, logl);
        //printf("qiter %d %f\n", i, logl);
        removeImpliedSpecNodes(tree, addedNodes);
        recon.setSize(tree->nnodes);
        events.setSize(tree->nnodes);
        
        Tree *tree2 = tree->copy();
        
        // save tree and logl
        trees.append(tree2);
        logls.append(logl);
        
        //if (isLogLevel(LOG_MEDIUM)) {
        //    printLog(LOG_MEDIUM, "qiter %d %f\n", i, logl);
        //    displayTree(tree2, getLogFile());
        //}
        
        proposer->revert(tree);
    }    
    
    // propose one of the subproposals 
    float choice = frand();
    float partsum = -INFINITY;
    
    for (int i=0; i<quickiter; i++) {
        partsum = logadd(partsum, logls[i]);
    
        //printf("part %d %f (%f)\n", i, expf(partsum - sum), choice);
        
        if (choice < expf(partsum - sum)) {
            // propose tree i
            // trees[i];
            //printLog(LOG_MEDIUM, "choose %d\n", i);
            tree->setTopology(trees[i]);
            
            //displayTree(trees[i]);
            //displayTree(tree);
            
            break;
        }
        
        
    }

    
    // clean up subproposals
    for (int i=0; i<quickiter; i++)
        delete trees[i];
    
    //proposer->propose(tree);
}

void DupLossProposer::revert(Tree *tree)
{
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->revert(tree);
        return;
    }
    
    //printf("set oldtop\n");
    tree->setTopology(oldtop);
    delete oldtop;
}


void DupLossProposer::reset()
{
    iter = 0;
}

void DupLossProposer::accept(bool accepted)
{
}



/*
bool treePropCmp(const DupLossProposer::TreeProp &a, 
                 const DupLossProposer::TreeProp &b)
{
    return a.second > b.second;
}


void DupLossProposer::queueTrees(Tree *tree)
{
    
    // recon tree to species tree
    recon.ensureSize(tree->nnodes);
    events.ensureSize(tree->nnodes);
    recon.setSize(tree->nnodes);
    events.setSize(tree->nnodes);
    
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    // init queue for subproposals
    clearQueue();
    //queue.clear();
    
    sum = -INFINITY;
    
    // make many subproposals
    proposer->reset();
    for (int i=0; i<quickiter; i++) {
        proposer->propose(tree);
        
        Node *oldroot1 = tree->root->children[0];
        Node *oldroot2 = tree->root->children[1];

        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        float logl = birthDeathTreePriorFull(tree, stree, recon, events, 
                                             dupprob, lossprob,
                                             doomtable, maxdoom);
        sum = logadd(sum, logl);
        
        // save tree and logl
        Tree *tree2 = tree->copy();
        queue.push_back(TreeProp(tree2, logl));

        // restore tree
        tree->reroot(oldroot1, oldroot2);
        proposer->revert(tree);
    }
    
    // makes the random sample slightly faster
    //std::sort(queue.begin(), queue.end(), treePropCmp);
    
    // start a new sample count
    samplei = 0;
}


void DupLossProposer::propose(Tree *tree)
{
    iter++;
    treesize = tree->nnodes;

    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->propose(tree);
        return;
    }
    
    if (queue.size() == 0 || samplei >= nsamples)
        queueTrees(tree);
    samplei++;
    
    // save old topology
    if (oldtop)
        delete oldtop;
    oldtop = tree->copy();
    
    // randomly sample a tree from the queue
    float choice = log(frand()) + sum;
    float partsum = -INFINITY;
    
    for (unsigned int i=0; i<queue.size(); i++) {
        partsum = logadd(partsum, queue[i].second);
        
        if (choice < partsum) {
            // propose tree i
            tree->setTopology(queue[i].first);

            // remove tree i from queue
            sum -= queue[i].second;
            queue[i] = queue.back();
            queue.pop_back();            
            break;
        }
    }
    
    testCorrect(tree); // debug: keep track of correct tree in search
}


void DupLossProposer::accept(bool accepted)
{
    //if (accepted) {
        clearQueue();
        
        // extend more iterations
        //if (extend)
        //    niter = max(niter, iter + treesize * 2);
        //}
}


void DupLossProposer::revert(Tree *tree)
{
    // do simple proposal if dup/loss probs are disabled
    if (dupprob < 0.0 || lossprob < 0.0 || quickiter <= 1) {
        proposer->revert(tree);
        return;
    }
    
    tree->setTopology(oldtop);
    delete oldtop;
    oldtop = NULL;
}

void DupLossProposer::reset()
{
    iter = 0;
    clearQueue();
}


void DupLossProposer::clearQueue()
{
    // clear queue
    while (queue.size() > 0) {
        Tree *tree = queue.back().first;
        queue.pop_back();
        delete tree;
    }
}
*/


//=============================================================================
// Fitting branch lengths


HkyFitter::HkyFitter(int nseqs, int seqlen, char **seqs, 
                     float *bgfreq, float tsvratio, int maxiter,
                     bool useLogl) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs),
    bgfreq(bgfreq),
    tsvratio(tsvratio),
    maxiter(maxiter),
    useLogl(useLogl)

{}

float HkyFitter::findLengths(Tree *tree)
{ 
    Timer timer;
    float logl = findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, 
                                        tsvratio, maxiter);
    runtime += timer.time();

    if (useLogl)
        return logl;
    else
        return 0.0;
}


//=============================================================================
// Prior function

SpidirPrior::SpidirPrior(
    int nnodes, SpeciesTree *stree, 
    SpidirParams *params, 
    int *gene2species,
    float predupprob, float dupprob, float lossprob, 
    int nsamples, bool approx, bool estGenerate) :
    
    nnodes(nnodes),
    stree(stree),
    params(params),
    gene2species(gene2species),
    recon(nnodes),
    events(nnodes),
    predupprob(predupprob),
    dupprob(dupprob),
    lossprob(lossprob),
    nsamples(nsamples),
    approx(approx),
    estGenerate(estGenerate)
{
    const int maxdoom = 20;

    doomtable = new float [stree->nnodes];
    calcDoomTable(stree, dupprob, lossprob, maxdoom, doomtable);
}


SpidirPrior::~SpidirPrior()
{
    delete [] doomtable;
}

float SpidirPrior::branchPrior(Tree *tree)
{
    Timer timer;

    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    float generate;
    
    if (estGenerate)
        generate = -1;
    else
        generate = -99;
        
    float logp = spidir::branchPrior(tree, stree,
                                     recon, events, params,
                                     generate, predupprob, dupprob, lossprob,
                                     nsamples, approx);
    branch_runtime += timer.time();

    return logp;
}



float SpidirPrior::topologyPrior(Tree *tree)
{
    Timer timer;

    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    
    const int maxdoom = 20;

    float logp = birthDeathTreePriorFull(tree, stree, recon, events, 
                                         dupprob, lossprob,
                                         doomtable, maxdoom);

    top_runtime += timer.time();

    return logp;
}


//=============================================================================


// propose initial tree by Neighbor Joining
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs)
{
    int nnodes = nseqs * 2 - 1;

    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);

    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    tree->setLeafNames(genes);

    return tree;
}


// propose initial tree and root by species
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species)
{
    Tree *tree = getInitialTree(genes, nseqs, seqlen, seqs);
    reconRoot(tree, stree, gene2species);

    return tree;
}


void printSearchStatus(Tree *tree, SpeciesTree *stree, int *gene2species,
                       int *recon=NULL, int *events=NULL)
{
    if (stree) {
        bool cleanupRecon = false;
        if (!recon) {
            recon = new int [tree->nnodes];
            cleanupRecon = true;
        }

        bool cleanupEvents = false;    
        if (!events) {
            events = new int [tree->nnodes];
            cleanupEvents = true;
        }    

        //assert(tree->assertTree());
    
        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        int losses = countLoss(tree, stree, recon);        
        
        // count dups
        int dups = 0;
        for (int i=0; i<tree->nnodes; i++)
            if (events[i] == EVENT_DUP)
                dups++;
        
        printLog(LOG_LOW, "search: dups = %d\n", dups);
        printLog(LOG_LOW, "search: loss = %d\n", losses);
        
        if (cleanupRecon)
            delete [] recon;
        if (cleanupEvents)
            delete [] events;        
    }                
    
    if (isLogLevel(LOG_LOW))
        displayTree(tree, getLogFile());
    
    if (isLogLevel(LOG_MEDIUM)) {
        printLog(LOG_MEDIUM, "tree: ");
        tree->writeNewick(getLogFile(), NULL, 0, true);
        printLog(LOG_MEDIUM, "\n\n");
    }
}


//=============================================================================
// Search Climb

TreeSearchClimb::TreeSearchClimb(Prior *prior,
				 TopologyProposer *proposer,
				 BranchLengthFitter *fitter) :
    prior(prior),
    proposer(proposer),
    fitter(fitter)
{
}


TreeSearchClimb::~TreeSearchClimb()
{}




Tree *TreeSearchClimb::search(Tree *initTree, 
			      string *genes, 
			      int nseqs, int seqlen, char **seqs)
{
    Tree *toptree = NULL;
    float toplogp = -INFINITY, nextlogp;
    Tree *tree = initTree;
    Timer correctTimer;
    
    // search testing
    Tree *correct = proposer->getCorrect();
    float correctLogp = -INFINITY;
    if (correct) {
        // determine probability of correct tree
        parsimony(correct, nseqs, seqs); // get initial branch lengths
        float seqlk = fitter->findLengths(correct);
        float branchp = prior->branchPrior(correct);
        float topp = prior->topologyPrior(correct);
        correctLogp = seqlk + branchp + topp;

        printLog(LOG_LOW, "search: correct tree lnl = %f\n", correctLogp);
    }

    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs,
                              prior->getSpeciesTree(), 
                              prior->getGene2species());
    
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // calc probability of initial tree
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    float seqlk = fitter->findLengths(tree);
    float branchp = prior->branchPrior(tree);
    float topp = prior->topologyPrior(tree);
    toplogp = seqlk + branchp + topp;
    toptree = tree->copy();
    
    
    // log initial tree
    printLog(LOG_LOW, "search: initial\n");
    printLog(LOG_LOW, "search: lnl    = %f\n", toplogp);
    printLog(LOG_LOW, "search: seqlk  = %f\n", seqlk);
    printLog(LOG_LOW, "search: branch = %f\n", branchp);          
    printLog(LOG_LOW, "search: top    = %f\n", topp);
    printSearchStatus(tree, prior->getSpeciesTree(), 
                      prior->getGene2species(), &recon[0], &events[0]);
    
    
    int accept = 0;
    int reject = 0;
    
    // search loop
    proposer->reset();
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposer->propose(tree);
        
	if (isLogLevel(LOG_MEDIUM)) {
	    displayTree(tree, getLogFile());
	}

        // calculate likelihood
        seqlk = fitter->findLengths(tree);
        branchp = prior->branchPrior(tree);
        topp = prior->topologyPrior(tree);
        nextlogp = seqlk + branchp + topp;
        

        // search test
        if (correct && (nextlogp >= correctLogp || 
                        tree->sameTopology(correct)))
        {
            printLog(LOG_LOW, "search: correct tree time = %f\n", 
                     correctTimer.time());
            printLog(LOG_LOW, "search: correct tree logl = %f\n", 
                     correctLogp);
            printLog(LOG_LOW, "search: correct tree is best = %d\n",
                     int(tree->sameTopology(correct)));
            printLog(LOG_LOW, "search: correct tree better logl = %f\n", 
                     nextlogp);

            correct = NULL;
        }

        // acceptance rule
        bool accept = (nextlogp > toplogp);

        // print accept
        if (accept)
            printLog(LOG_LOW, "search: accept\n");
        else
            printLog(LOG_LOW, "search: reject\n");


        // print info
        printLog(LOG_LOW, "search: lnl    = %f\n", nextlogp);
        printLog(LOG_LOW, "search: seqlk  = %f\n", seqlk);
        printLog(LOG_LOW, "search: branch = %f\n", branchp);          
        printLog(LOG_LOW, "search: top   = %f\n", topp);
        

        // act on acceptance
        if (accept) {
            // accept
            accept++;
            proposer->accept(true);
            toplogp = nextlogp;
            delete toptree;
            toptree = tree->copy();

            printSearchStatus(tree, 
                              prior->getSpeciesTree(), 
                              prior->getGene2species(), &recon[0], &events[0]);

        } else {           
            // reject, undo topology change
            
            // display rejected tree
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
				  prior->getGene2species(), 
                                  &recon[0], &events[0]);
            
            reject++;
            proposer->accept(false);            
            proposer->revert(tree);
        }
    }
    
    printLog(LOG_LOW, "accept rate: %f\n", accept / float(accept+reject));
    
    return toptree;
}




extern "C" {


// Calculate the likelihood of a tree
Tree *searchClimb(int niter, int quickiter,
		  int nseqs, char **gene_names, char **seqs,
		  int nsnodes, int *pstree, float *sdists,
		  int *gene2species,
		  float *sp_alpha, float *sp_beta, float generate,
		  float pretime_lambda, float birth, float death,
		  float gene_alpha, float gene_beta,
		  float *bgfreq, float kappa,
		  int nsamples, bool approx)
{
    int nnodes = 2*nseqs - 1;

    setLogLevel(LOG_MEDIUM);

    // create stree
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    stree.setDists(sdists);
    
    // params
    SpidirParams params(nsnodes, NULL, sp_alpha, sp_beta, 
			gene_alpha, gene_beta, pretime_lambda);
    
    // make gene names
    string genes[nseqs];
    //char s[101];
    for (int i=0; i<nseqs; i++) {
	//snprintf(s, 100, "%d", i);
	genes[i] = gene_names[i];
    }

    // prior
    SpidirPrior prior(nnodes, &stree, &params, gene2species,
		      pretime_lambda, birth, death, nsamples, approx, false);

    // seq likelihood
    const int maxiter = 10;
    int seqlen = strlen(seqs[0]);
    HkyFitter fitter(nseqs, seqlen, seqs, 
		     bgfreq, kappa, maxiter);

    // proposers
    SprNbrProposer proposer2(&stree, gene2species, niter);
    DupLossProposer proposer(&proposer2, &stree, gene2species, 
			     birth, death,
                             quickiter, niter);
    
    TreeSearchClimb search(&prior, &proposer, &fitter);

    Tree *tree = search.search(NULL, genes, nseqs, seqlen, seqs);
    
    return tree;
}

} // extern "C"






//=============================================================================
// MCMC search


Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 SampleFunc *samples,
                 Prior *prior,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter)
{
    Tree *toptree = NULL;
    float toplogl = -INFINITY, logl=-INFINITY, nextlogl;
    Tree *tree = initTree;
    
    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs);
    
    // used for printSearchStatus
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // init likelihood score
    parsimony(tree, nseqs, seqs); // get initial branch lengths
    logl = fitter->findLengths(tree);
    logl += prior->branchPrior(tree);
    toplogl = logl;
    toptree = tree->copy();
    
    int accept = 0;
    int reject = 0;
    
        
    // MCMC loop
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposer->propose(tree);
        //assert(tree->assertTree());
        
        float seqlk = fitter->findLengths(tree);
        float branchlk = prior->branchPrior(tree);
        nextlogl = seqlk + branchlk;       
        
        
        // acceptance rule
        if (nextlogl > logl ||
            nextlogl - logl > log(frand()))
        {
            // accept
            printLog(LOG_MEDIUM, "search: accept %f (%f)\n", nextlogl, logl);
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
                                 prior->getGene2species(), recon, events);            
            
            accept++;
            logl = nextlogl;

            // keep track of toptree            
            if (logl > toplogl) {
                delete toptree;
                toptree = tree->copy();
                toplogl = logl;
            }
            
        } else {                
            // reject, undo topology change
            printLog(LOG_MEDIUM, "search: reject %f (%f)\n", nextlogl, logl);
            if (isLogLevel(LOG_MEDIUM))
                printSearchStatus(tree, prior->getSpeciesTree(), 
                                  prior->getGene2species(), recon, events);
            
            reject++;
            proposer->revert(tree);
        }
        
        // return a tree sample
        (*samples)(tree);
    }
    
    printLog(LOG_LOW, "accept rate: %f\n", accept / float(accept+reject));
    
    return toptree;
}



} // namespace spidir
