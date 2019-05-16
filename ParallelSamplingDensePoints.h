//NUEVO REPO
#include <stdio.h>
#include <time.h>
#include <limits>
#include <vector>
#include <list>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include "sample.h"

#define N_THREADS 2
#define WIDTH 500
#define HEIGHT 500
#define L 4000

using namespace std;


int numberPoints;
int r = 30;
float size_ = r / sqrt(2);
int cols = WIDTH / size_;
int rows = HEIGHT / size_;//rows

vector<Sample *> cloud(L); //dense points cloud
vector<int> pointIndex(L);
vector<list<Sample *>> grid(cols *rows); //inicializar el grid !
omp_lock_t writelock;



class ParallelSamplingDensePoints
{
    Sample **S = new Sample *[cols * rows];

   
    
public:
    

    void printGrid()
    {   
        std::cout<<L<<'\n';
        for (int i =0;i<rows*cols;i++)
        {
            if (S[i]->status!="0")
            {
                cout<<S[i]->pos[0]<<" "<<S[i]->pos[1]<<endl;
            }
        }
    }

    bool existeIddle()
    {
        for (auto &p : cloud)
        {
            if (p->status=="IDDLE")
            {
                return true;
            }
        }
        return false;
    }

    void detectCollision(Sample *p, int r)
    {
        int col = ((p->pos[0]) / size_);
        int row = ((p->pos[1]) / size_);

        

        for (int i = -(r-1); i <= (r-1); i++)
        {
            for (int j = -(r-1); j <= (r-1); j++)
            {
                int position = ((col + i) + (row + j)) * cols;
                if (position >= 0 && position < (cols * rows))
                {
                    std::list<Sample *> neighbour = grid[position];
                    for (auto &sample : neighbour)
                    {
                        if (sample->status == "IDDLE")
                        {
                            p->I.push_back(sample);
                        }
                        else if (sample->status == "ACTIVE")
                        {
                            p->A.push_back(sample);
                        }
                    }
                }
            }
        }
        int sad = (col + row) * rows;
        for (auto &samples : grid[sad])
        {
            if (samples->pos[0] != p->pos[0] && samples->pos[1] != p->pos[1])
            {
                if (samples->status == "IDDLE")
                {
                    p->I.push_back(samples);
                }
                else if (samples->status == "ACTIVE")
                {
                    p->A.push_back(samples);
                }
            }
        }
    }
    void checkStatus(Sample *pi)//pi as argument
    {
        omp_set_lock(&writelock);
        if(pi->status != "ACTIVE"){//atomic
            omp_unset_lock(&writelock);
            return ;
        }
        omp_unset_lock(&writelock);

        for (auto &q :pi->A) {
            if(q->prioridad > pi->prioridad){
                checkStatus(q);
                
                omp_set_lock(&writelock);
                if (q->status=="ACCEPTED"){//atomic
                    pi->status = "REJECTED";
                    omp_unset_lock(&writelock);
                    return;
                }
                omp_unset_lock(&writelock);
            }
        }
#pragma omp critical
        {
            pi->status = "ACCEPTED";//atomic
        }
        return;
    }

    void generateDensePoints()
    {
        //printf("puntos en la nube : %lu\n",cloud.size());
        random_shuffle(pointIndex.begin(), pointIndex.end());
        //int count = 0;
        //std::cout<<L<<'\n';
        //for(auto i : grid){
        //    for(auto j : i){
        //        std::cout<<j->pos[0]<<" "<<j->pos[1]<<std::endl;
        //    }
        //}
    }
    
    void parallel_generateDP()
    {



        omp_set_num_threads(N_THREADS);
        omp_set_num_threads(N_THREADS);
#pragma omp parallel for
            for(int i = 0; i < rows*cols; ++i)
            {
               Sample *sample = new Sample();
               sample->status = "0";
               S[i]=sample;
            }
            omp_set_num_threads(N_THREADS);
#pragma omp parallel for
            for (int i = 0; i < L; ++i)
            {
                int id = omp_get_thread_num();
                //printf("%d\n", id);
                Sample *sample = new Sample();
                sample->pos[0] = rand() % (WIDTH - 1); //numeros aleatorios entre 0 y 500
                sample->pos[1] = rand() % (HEIGHT - 1);
                sample->status = "IDDLE";
                //  printf("x : %d y : %d\n", sample->pos[0], sample->pos[1]);
                cloud[i] = sample;
                //agregar al grid cell correspondiente
                int ii = sample->pos[0] / size_;
                int j = sample->pos[1] / size_;
                int id_grid_cell = ii + j * cols;
                sample->id_grid_cell = id_grid_cell;
                // printf("id : %d cell: %d \n",data1->id_t, id_grid_cell);
                grid[id_grid_cell].push_back(sample);
                pointIndex[i] = i;
            }

    };

    void run_parallel()
    {
        int i =0;
        do {
            omp_set_num_threads(N_THREADS);
#pragma omp parallel num_threads(N_THREADS)
            {
                int id_thread = omp_get_thread_num() + 1;
                int start = (id_thread - 1) / N_THREADS * L;
                int fin = (id_thread / N_THREADS) * L;
                int curr_index = start;

                while (cloud[pointIndex[curr_index]]->status != "IDDLE" ) {
                    curr_index++;
                    if (curr_index == fin) break;
                }
                Sample *sample = cloud[pointIndex[curr_index]];
                sample->status = "ACTIVE";
                sample->prioridad = (rand() * N_THREADS + id_thread) / numeric_limits<int>::max();
#pragma omp barrier
                detectCollision(sample, 2 * r);
#pragma omp barrier
                checkStatus(sample);
                if (sample->status == "ACCEPTED") {
                    omp_set_lock(&writelock);
                    S[sample->id_grid_cell] = sample; //atomic
                    omp_unset_lock(&writelock);

                    for (auto &i :sample->I) {
                        omp_set_lock(&writelock);
                        i->status = "REJECTED";
                        omp_unset_lock(&writelock);
                    }
                }
            }
            i++;
       // }while(i<500);
        }while(existeIddle());
        
    }

    
    
};
