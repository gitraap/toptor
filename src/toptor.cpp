#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
//#include <gtest/gtest.h>

//#define TOPGEN_PBGL //manual definition

#include <topgen.hpp>

#ifdef TOPGEN_PBGL
#include <boost/graph/use_mpi.hpp>
#endif 

//Time
#include <omp.h>
#include <sys/time.h>

// choose only one
//#define REAL_TIMER
#define VIRTUAL_TIMER
//#define PROF_TIMER
#define MYTIME
#ifdef MYTIME
#define DECLARE(vble) struct timespec vble
#define DELAY(vble,begin,end) ( vble = (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / 1000000000.0 )
#define INCDELAY(vble,begin,end) ( vble += (end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec) / 1000000000.0 )
#define READ_CLOCK(c) (clock_gettime(CLOCK_MONOTONIC, &c))
#else
#define DECLARE(vble) ;
#define DELAY(vble,begin,end) ;
#define INCDELAY(vble,begin,end) ;
#define READ_CLOCK(c) ;
#endif

unsigned int getTotalAllocMemory() {
    //Note: this value is in KB!
    FILE * file = fopen("/proc/self/status", "r");
    if (file == NULL) {
        perror("/proc/self/status");
    }

    unsigned int result = 0;
    char line[128];

    while (fgets(line, 128, file) != NULL) {
        //printf("LINE:: %s\n", line);
        if (1 == sscanf(line, "VmPeak: %u ", &result)) {
            //printf("%s ==> %u\n", line, result);
            break;
        }
    }

    // HACK /proc/self/status cannot be close
    /*if (fclose(file) != 0)
        perror("/proc/self/status");*/

    return result;
}

//#define PRINT_TOPOLOGYMANAGER_INFORMATION

class CParameters {
public:
    /* Topology Generator */
    topgen::topologyType m_eTopology;
    unsigned int m_iEndNodes;
    unsigned int m_iRadix;
    int * m_pRadixes;
    unsigned int m_iInternalPorts;
    topgen::connectionType m_eConnection;
    topgen::rtAlgorithmType m_eRoutingAlgorithm;
    topgen::nodeType m_eEmbedded;
    topgen::matchingType m_eMatching;
    unsigned int m_iKarity1;
    unsigned int m_iKarity2;
    unsigned int m_iKarity3;
    unsigned int m_iDimensions;
    int * m_pKaries;
    int m_iBits;
    char * m_sPatternFile;
    int m_iHeight;
    int * m_pChildren;
    int * m_pParents;
    topgen::KNSIndirectEnum m_tIndirectTopo;
    bool m_bVerbose;
    int m_iTopologyParams;
    char ** m_pTopologyParams;
    int m_iRoutingParams;
    char ** m_pRoutingParams;
    int m_iCounter;
    int dfly_a;
    int dfly_h;
    int dfly_p;
    bool m_bDotFlag;
    char * m_pDotFilename;
    
    CParameters() {
        m_pRadixes = NULL;
        m_pChildren = NULL;
        m_pParents = NULL;
        m_iCounter = 1; // the command (i.e. ./toptor)
        /* DO NOT MODIFY THIS VALUES */
        m_iInternalPorts = 0;
        m_eConnection = topgen::butterfly;
        m_eRoutingAlgorithm = topgen::destro;
        m_eEmbedded = topgen::switchSTD;
        m_eMatching = topgen::Alpha;
        m_pKaries = NULL;
        m_sPatternFile = NULL;
        m_bVerbose = false;
        m_bDotFlag = false;
        m_pDotFilename = NULL;
    };

    ~CParameters() {
        if (m_pKaries)
            delete [] m_pKaries;
        if (m_sPatternFile)
            delete [] m_sPatternFile;
        if (m_pRadixes)
            delete [] m_pRadixes;
    };
};

// Global variables
int gHEIGHT;
int gNICs_MIN;
int gNICs_MAX;
int gRADIX;

void solveparents(int productorio, int * m, int j, int * w) {
    if (j == gHEIGHT + 1) {
        bool solution_found = true;
        // check the sanity of this configuration
        if (m[gHEIGHT] == gRADIX && w[1] == 1) {
            // the components are not null
            for (int jj = 1; jj <= gHEIGHT; jj++) {
                if (m[jj] == 0 || w[jj] == 0) {
                    solution_found = false;
                    break;
                }
            }
            // the addition of any switch port at a given level is equal to the radix (level 0)
            for (int jj = 1; jj + 1 <= gHEIGHT; jj++) {
                if (((m[jj] + w[jj + 1]) == gRADIX) == false)
                    solution_found = false;
            }
        } else {
            solution_found = false;
        }

        if (solution_found) {
            printf("N %d :: XGFT (h=%d; ", productorio, gHEIGHT);
            for (int jj = 1; jj <= gHEIGHT; jj++) {
                if (jj > 1) printf(",");
                printf("%d", m[jj]);
            }
            printf("; ");
            for (int jj = 1; jj <= gHEIGHT; jj++) {
                if (jj > 1) printf(",");
                printf("%d", w[jj]);
            }
            printf(")");

            // get information of the type of XGFT
            bool is_slim = true;
            bool is_fat = true;
            bool is_full = true;
            for (int jj = 1; jj + 1 <= gHEIGHT; jj++) {
                if ((m[jj] >= w[jj + 1]) == false)
                    is_slim = false;
                if ((m[jj] <= w[jj + 1]) == false)
                    is_fat = false;
                if ((m[jj] == w[jj + 1]) == false)
                    is_full = false;
            }
            if (is_slim) printf(" Slim ");
            else printf(" ---- ");
            if (is_fat) printf(" Fat  ");
            else printf(" ---- ");
            if (is_full) printf(" Full ");
            else printf(" ---- ");

            // print ratios between 1:1 and 2:1 (Xin Yuan states they are the best)
            // calculate standard deviation to identify the best ratios
            bool good_ratio = true;
            double sum = 0.0;
            double * x = new double[gHEIGHT];
            for (int jj = 1; jj + 1 <= gHEIGHT; jj++) {
                double ratio = (double) m[jj] / w[jj + 1];
                sum += ratio;
                x[jj - 1] = ratio;

                printf("RT%d %.2f ", jj, ratio);
                if ((1.0 <= ratio && ratio <= 2.0) == false) {
                    good_ratio = false;
                    break;
                }
            }

            int n = gHEIGHT;
            double mean = sum / n;
            double dev, sdev = 0.0;

            for (int i = 0; i < n; ++i) {
                dev = (x[i] - mean)*(x[i] - mean);
                sdev = sdev + dev;
            }
            delete [] x;

            double var = sdev / (n - 1);
            double sd = sqrt(var);
            double cv = (sd / mean) * 100;

            if (good_ratio) printf(" GoodRatio! VAR %.3f SD %.3f CV %.3f", var, sd, cv);

            printf("\n");

        }
    } else {
        /*for (int r = 0; r <= gRADIX; r++) {
        //for (int r = 1; r <= gRADIX - m[j]; r++) {
            int * new_w = new int[gHEIGHT + 1];
            for (int jj = 1; jj <= gHEIGHT; jj++)
                new_w[jj] = w[jj];
            new_w[j] = r;

            solveparents(productorio, m, j + 1, new_w);
            delete [] new_w;
        }*/
        int * new_w = new int[gHEIGHT + 1];
        for (int jj = 1; jj <= gHEIGHT; jj++)
            new_w[jj] = w[jj];

        new_w[j] = (j == 1) ? 1 : gRADIX - m[j - 1];

        solveparents(productorio, m, j + 1, new_w);
        delete [] new_w;
    }
}

void solvechildren(int j, int * m) {
    if (j == gHEIGHT + 1) {
        int productorio = 1;
        for (int jj = 1; jj <= gHEIGHT; jj++) {
            productorio *= m[jj];
        }

        /*printf("productorio %d W[", productorio);
        for (int jj = 1; jj <= gHEIGHT; jj++) {
            if (jj > 1) printf(",");
            printf("%d", m[jj]);
        }
        printf("]\n");*/

        if (gNICs_MIN <= productorio && productorio <= gNICs_MAX) {
            int * Ws = new int[gHEIGHT + 1]; // parents. omit position 0
            for (int i = 1; i <= gHEIGHT; i++) {
                Ws[i] = 0;
            }

            solveparents(productorio, m, 1, Ws);

            delete [] Ws;
        }
    } else {
        for (int r = 0; r <= gRADIX; r++) {
            int * new_m = new int[gHEIGHT + 1];
            for (int jj = 1; jj <= gHEIGHT; jj++)
                new_m[jj] = m[jj];
            new_m[j] = r;

            solvechildren(j + 1, new_m);
            delete [] new_m;
        }
    }
}

int main(int argc, char* argv[]) {
    
    if (argc == 1) {       
#ifndef TOPGEN_PBGL   
        printf("usage: [best] topology [options] [memory] [bgl options] [memory]\n");
        printf("\n");
        printf("Available options for topology:\n");
        printf("  fully nodes\n");
        printf("  router nodes\n");
        printf("  custom1\n");
        printf("  custom2\n");
        printf("  cube dimensions\n");
        printf("  min endnodes radix\n");
        printf("  mesh dimensions k1 k2 k3\n");
        printf("  torus dimensions k1 k2 k3\n");
        printf("  torusnd dimensions k1 k2 ... kn\n");
        printf("  kns dimensions k1 k2 ... kN router|min radix1 radix2 ... radixN [verbose]\n");
        printf("  xgft height M1 ... Mh W1 ... Wh radix [verbose]\n");
        printf("  rlft end-nodes radix [verbose]\n");
        printf("  dfly a h p [verbose]\n");
        printf("\n");
        printf("  best xgft min_endnodes max_endnodes radix max_height\n");
        printf("\n");
        printf("Available options for bgl\n");
        printf("  graphviz outfile.gv\n");
        printf("  graphml  outfile.gm\n");
        printf("  bfs from\n");
        printf("  dfs from\n");
        printf("  dijkstra from to\n");
        printf("  bellman-ford from to\n");
        printf("  prim from\n");
        printf("  kruskal from\n");
#else
        printf("usage: [memory] [pbgl options] [memory]\n");
        printf("\n");
        printf("Available options for pbgl\n");
        printf("  graphviz infile.gv\n");
        printf("  dijkstra from to\n");
#endif
        
        printf("\n");
        return 0;
    }
    
    //omp_set_num_threads(4);
    //std::cout << "+++ omp threads=" << omp_get_num_threads() << std::endl;

    /* Initialization of the library */
    CParameters * g_pParameter = new CParameters();
    topgen::topology * m_pTopology;
    
#ifndef TOPGEN_PBGL
    if (!strcmp(argv[1], "fully")) {
        g_pParameter->m_eTopology = topgen::fully;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iEndNodes = atoi(argv[g_pParameter->m_iCounter++]);
        m_pTopology = new topgen::CFully(g_pParameter->m_iEndNodes);
    } else if (!strcmp(argv[1], "router")) {
        g_pParameter->m_eTopology = topgen::router;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iEndNodes = atoi(argv[g_pParameter->m_iCounter++]);
        m_pTopology = new topgen::CRouter(g_pParameter->m_iEndNodes);
    } else if (!strcmp(argv[1], "custom1")) {
        g_pParameter->m_eTopology = topgen::custom1;
        g_pParameter->m_iCounter++;
        m_pTopology = new topgen::CCustom1();
    } else if (!strcmp(argv[1], "custom2")) {
        g_pParameter->m_eTopology = topgen::custom2;
        g_pParameter->m_iCounter++;
        m_pTopology = new topgen::CCustom2();
    } else if (!strcmp(argv[1], "mesh")) {
        g_pParameter->m_eTopology = topgen::mesh;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iDimensions = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iKarity1 = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iKarity2 = atoi(argv[g_pParameter->m_iCounter++]);
        if (g_pParameter->m_iDimensions == 3) {
            g_pParameter->m_iKarity3 = atoi(argv[g_pParameter->m_iCounter++]);
        }
        m_pTopology = new topgen::CMesh(g_pParameter->m_iDimensions, g_pParameter->m_iKarity1, g_pParameter->m_iKarity2, g_pParameter->m_iKarity3);
    } else if (!strcmp(argv[1], "torus")) {
        g_pParameter->m_eTopology = topgen::torus;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iDimensions = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iKarity1 = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iKarity2 = atoi(argv[g_pParameter->m_iCounter++]);
        if (g_pParameter->m_iDimensions == 3) {
            g_pParameter->m_iKarity3 = atoi(argv[g_pParameter->m_iCounter++]);
        }
        m_pTopology = new topgen::CTorus(g_pParameter->m_iDimensions, g_pParameter->m_iKarity1, g_pParameter->m_iKarity2, g_pParameter->m_iKarity3);
        ((topgen::CTorus*)m_pTopology)->Settings(g_pParameter->m_eEmbedded);
    } else if (!strcmp(argv[1], "cube")) {
        g_pParameter->m_eTopology = topgen::cube;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iDimensions = atoi(argv[g_pParameter->m_iCounter++]);
        m_pTopology = new topgen::CCube(g_pParameter->m_iDimensions);
    } else if (!strcmp(argv[1], "min")) {
        g_pParameter->m_eTopology = topgen::min;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iEndNodes = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iRadix = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_eConnection = topgen::butterfly;
        g_pParameter->m_eRoutingAlgorithm = topgen::destro;
        m_pTopology = new topgen::CMIN(g_pParameter->m_iEndNodes, g_pParameter->m_iRadix, g_pParameter->m_eConnection, g_pParameter->m_eRoutingAlgorithm);
        ((topgen::CMIN*)m_pTopology)->Settings(topgen::standard, g_pParameter->m_eEmbedded, g_pParameter->m_eMatching, g_pParameter->m_iInternalPorts, g_pParameter->m_sPatternFile);
    } else if (!strcmp(argv[1], "torusnd")) {
        g_pParameter->m_eTopology = topgen::torusND;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iDimensions = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_pKaries = new int[g_pParameter->m_iDimensions];
        for (unsigned int d = 0; d < g_pParameter->m_iDimensions; d++) {
            g_pParameter->m_pKaries[d] = atoi(argv[g_pParameter->m_iCounter++]);
        }
        m_pTopology = new topgen::CTorusND(g_pParameter->m_iDimensions, g_pParameter->m_pKaries);
    } else if (!strcmp(argv[1], "kns")) {
        g_pParameter->m_eTopology = topgen::kns;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iDimensions = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_pKaries = new int[g_pParameter->m_iDimensions];
        unsigned int dim;
        for (unsigned dim = 0; dim < g_pParameter->m_iDimensions; dim++) {
            g_pParameter->m_pKaries[dim] = atoi(argv[g_pParameter->m_iCounter++]);
        }

        if (!strcmp(argv[g_pParameter->m_iCounter], "min")) {
            g_pParameter->m_tIndirectTopo = topgen::knsmin;
            g_pParameter->m_iCounter++;
        } else if (!strcmp(argv[g_pParameter->m_iCounter], "router")) {
            g_pParameter->m_tIndirectTopo = topgen::knsrouter;
            g_pParameter->m_iCounter++;
        } else {
            printf("Indirect Topology '%s' not supported", argv[g_pParameter->m_iCounter]);
            return -1;
        }

        g_pParameter->m_pRadixes = new int[g_pParameter->m_iDimensions];
        for (dim = 0; dim < g_pParameter->m_iDimensions; dim++) {
            g_pParameter->m_pRadixes[dim] = atoi(argv[g_pParameter->m_iCounter++]);
        }

        if (argv[g_pParameter->m_iCounter] != NULL &&
            !strcmp(argv[g_pParameter->m_iCounter], "verbose")) {
            g_pParameter->m_bVerbose = true;
            g_pParameter->m_iCounter++;
        } else {
            g_pParameter->m_bVerbose = false;
        }

        m_pTopology = new topgen::CKNS(g_pParameter->m_iDimensions,
            g_pParameter->m_pKaries,
            g_pParameter->m_tIndirectTopo,
            g_pParameter->m_pRadixes,
            g_pParameter->m_bVerbose);
    } else if (!strcmp(argv[1], "xgft")) {
        g_pParameter->m_eTopology = topgen::xgft;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iHeight = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_pChildren = new int[g_pParameter->m_iHeight];
        for (int children = 0; children < g_pParameter->m_iHeight; children++) {
            g_pParameter->m_pChildren[children] = atoi(argv[g_pParameter->m_iCounter++]);
        }
        g_pParameter->m_pParents = new int[g_pParameter->m_iHeight];
        for (int parents = 0; parents < g_pParameter->m_iHeight; parents++) {
            g_pParameter->m_pParents[parents] = atoi(argv[g_pParameter->m_iCounter++]);
        }

        g_pParameter->m_iRadix = atoi(argv[g_pParameter->m_iCounter++]);
        if (argv[g_pParameter->m_iCounter] != NULL &&
            !strcmp(argv[g_pParameter->m_iCounter], "verbose")) {
            g_pParameter->m_bVerbose = true;
            g_pParameter->m_iCounter++;
        } else {
            g_pParameter->m_bVerbose = false;
        }

        m_pTopology = new topgen::CXGFT(g_pParameter->m_iRadix,
            g_pParameter->m_iHeight,
            g_pParameter->m_pChildren,
            g_pParameter->m_pParents,
            g_pParameter->m_bVerbose);
    } else if (!strcmp(argv[1], "rlft")) {
        g_pParameter->m_eTopology = topgen::rlft;
        g_pParameter->m_iCounter++;
        g_pParameter->m_iEndNodes = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->m_iRadix = atoi(argv[g_pParameter->m_iCounter++]);
        if (argv[g_pParameter->m_iCounter] != NULL &&
            !strcmp(argv[g_pParameter->m_iCounter], "verbose")) {
            g_pParameter->m_bVerbose = true;
            g_pParameter->m_iCounter++;
        } else {
            g_pParameter->m_bVerbose = false;
        }

        m_pTopology = new topgen::CRLFT(g_pParameter->m_iEndNodes,
            g_pParameter->m_iRadix,
            topgen::butterfly,
            topgen::destro,
            g_pParameter->m_bVerbose);
    } else if (!strcmp(argv[1], "dfly")) {
        g_pParameter->m_eTopology = topgen::dfly;
        g_pParameter->m_iCounter++;
        g_pParameter->dfly_a = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->dfly_h = atoi(argv[g_pParameter->m_iCounter++]);
        g_pParameter->dfly_p = atoi(argv[g_pParameter->m_iCounter++]);
        if (argv[g_pParameter->m_iCounter] != NULL &&
            !strcmp(argv[g_pParameter->m_iCounter], "verbose")) {
            g_pParameter->m_bVerbose = true;
            g_pParameter->m_iCounter++;
        } else {
            g_pParameter->m_bVerbose = false;
        }

        m_pTopology = new topgen::CDFly(g_pParameter->dfly_a,
            g_pParameter->dfly_h,
            g_pParameter->dfly_p,
            topgen::dfly_canonic,
            topgen::minimal);
    } else if (!strcmp(argv[1], "best")) {
        g_pParameter->m_iCounter++;
        if (!strcmp(argv[g_pParameter->m_iCounter], "xgft")) {
            g_pParameter->m_iCounter++;
            gNICs_MIN = atoi(argv[g_pParameter->m_iCounter++]);
            gNICs_MAX = atoi(argv[g_pParameter->m_iCounter++]);
            gRADIX = atoi(argv[g_pParameter->m_iCounter++]);
            int maxHeight = atoi(argv[g_pParameter->m_iCounter++]);
            //const int maxRadix = 64;

            for (gHEIGHT = 2; gHEIGHT <= maxHeight; gHEIGHT++) {
                printf("***Height %d***\n", gHEIGHT);
                int * Ms = new int[gHEIGHT + 1]; // children. omit position 0
                for (int i = 0; i <= gHEIGHT; i++) {
                    Ms[i] = 1;
                }

                solvechildren(1, Ms);

                delete [] Ms;
            }

        } else {
            printf("Topology not supported");
            delete g_pParameter;
            return -1;
        }
        return 0;
    
    } else {
        printf("Topology not supported");
        delete g_pParameter;
        return -1;
    }

    double m_time_runall = 0.0;
    DECLARE(m_time_runall_b);
    DECLARE(m_time_runall_e);
    READ_CLOCK(m_time_runall_b);
    m_pTopology->RunAll();
    READ_CLOCK(m_time_runall_e);
    INCDELAY(m_time_runall, m_time_runall_b, m_time_runall_e); 
    
    std::cout << "+++ RunAll Time " << m_time_runall << " seconds" << std::endl;

    if (g_pParameter->m_bVerbose) {
        //printf("toptor :=> Printing your final topology ;)\n");
        //m_pTopology->PrintNetwork(stdout);
        m_pTopology->printTopology();

        printf("toptor :=> Printing your final routes ;)\n");
        //m_pTopology->PrintRoutes(stdout); ??
        
    }

    // libraries shouldn't print to output standard
    //printf("toptor :=> Printing your final paths ;)\n");
    //m_pTopology->PrintPaths(stdout);

    //topgen::g_pLog->Dump(stdout);
    
#endif

    if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "memory")) {
        g_pParameter->m_iCounter++;
        //m_pTopology->Outline();
        printf("+++ VmPeak1 (kB): %u\n", getTotalAllocMemory());
    }
    
#ifndef TOPGEN_PBGL
    
    if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "bgl")) {
        printf("*** BGL ***\n");
        g_pParameter->m_iCounter++;
        m_pTopology->build_bgl();
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "graphviz")) {
            printf("*** graphviz ***\n");
            g_pParameter->m_iCounter++;
            std::string outfile(argv[g_pParameter->m_iCounter++]);
            m_pTopology->write_graphviz(outfile);
            m_pTopology->run_read_graphviz(outfile); //it creates bgl-compare.gv
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "graphml")) {
            printf("*** graphml ***\n");
            g_pParameter->m_iCounter++;
            std::string outfile(argv[g_pParameter->m_iCounter++]);
            m_pTopology->write_graphml(outfile);
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "bfs")) {
            printf("*** Breadth-first search ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_breadth_first_search(from);
        }

        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "dfs")) {
            printf("*** Depth-first search ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_depth_first_search(from);
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "dijkstra")) {
            printf("*** dijkstra ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            const int to = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_dijkstra_shortest_paths(from, to);
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "bellman-ford")) {
            printf("*** bellman-ford ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            const int to = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_bellman_ford_shortest_paths(from, to);
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "prim")) {
            printf("*** prim ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_prim_minimum_spanning_tree(from);
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "kruskal")) {
            printf("*** kruskal ***\n");
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            m_pTopology->run_kruskal_minimum_spanning_tree(from);
        }
    }
        
#else //TOPGEN_PBGL
    
    printf("*** Loading PBGL environment ***\n");
    boost::mpi::environment env(argc, argv);
    
    boost::mpi::communicator world;
        
    std::cout << "*** Hello I am the rank " << world.rank() << " of " << world.size() << " ***\n";
    
    if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "pbgl")) {
        printf("*** PBGL ***\n");
        g_pParameter->m_iCounter++;
        std::string infile;
        // TODO improve this. Now create a network to have an object
        g_pParameter->m_eTopology = topgen::pbgl;
        m_pTopology = new topgen::CPBGL();
        printf("- Topology created\n");
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "graphviz")) {
            printf("-- option: graphviz begins\n");
            g_pParameter->m_iCounter++;
            std::string gvfile(argv[g_pParameter->m_iCounter++]);
            infile = gvfile;
            std::cout << "-- graphviz ends infile=" << infile << "\n";
        }
        
        if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "dijkstra")) {
            g_pParameter->m_iCounter++;
            const int from = atoi(argv[g_pParameter->m_iCounter++]);
            const int to = atoi(argv[g_pParameter->m_iCounter++]);
            printf("--- option: dijkstra from %d to %d ***\n",from,to);
            std::cout << "--- run_pbgl begins " << infile << ":" << from << ":" << to << "\n";
            ((topgen::CPBGL*)m_pTopology)->run_pbgl(infile, from, to);
            std::cout << "--- run_pbgl ends   " << infile << ":" << from << ":" << to << "\n";
        }
    }
    
#endif
    
    if (g_pParameter->m_iCounter < argc && !strcmp(argv[g_pParameter->m_iCounter], "memory")) {
        g_pParameter->m_iCounter++;
        //m_pTopology->Outline();
        printf("+++ VmPeak2 (kB): %u\n", getTotalAllocMemory());
    }
        
    delete m_pTopology;
    delete g_pParameter;
    return 0;
}
