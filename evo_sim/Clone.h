//
//  Clone.hpp
//  variable fitness branching process
//
//  Created by Debra Van Egeren on 2/17/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#ifndef clone_h
#define clone_h

#include <stdio.h>
#include <vector>
#include <queue>
#include <string>

using namespace std;

class CList;
class CellType;
class MutationHandler;

class Clone{
private:
    Clone *next_node;
    Clone *prev_node;
protected:
    long long cell_count;
    CellType *cell_type;
    
    double birth_rate;
    double mut_prob;
    
    
    virtual bool checkRep(){
        return !(mut_prob < 0 || birth_rate < 0 || cell_count < 0);
    }
public:
    /* adjusts the linked list so the clone is no longer in the list. adjusts root and end of cell_type as appropriate.
     MODIFIES cell_type
     */
    virtual ~Clone();
    
    Clone(CellType& type);
    
    Clone(CellType& type, double mut);
    
    virtual void reproduce() = 0;
    
    virtual void update(double t){}
    
    virtual bool hasDied(){
        return false;
    }
    
    virtual bool hasReproduced(){
        return false;
    }
    
    virtual bool readLine(vector<string>& parsed_line) = 0;
    
    void addCells(int num_cells);
    
    double getMutProb(){
        return mut_prob;
    }
    
    double getBirthRate(){
        return birth_rate;
    }
    double getDeathRate();
    double getTotalBirth(){
        return birth_rate * cell_count;
    }
    bool isSingleCell(){
        return cell_count == 1;
    }
    long long getCellCount(){
        return cell_count;
    }
    void setNext(Clone *child){
        next_node = child;
    }
    void setPrev(Clone *parent){
        prev_node = parent;
    }
    Clone& getNextWithinType(){
        return *next_node;
    }
    Clone& getPrevWithinType(){
        return *prev_node;
    }
    CellType& getType(){
        return *cell_type;
    }
    
    Clone* getNextClone();
    
    /* removes one cell from this clone's population
     should not be called if there is <=1 cell left in the clone
     MODIFIES clone_list, cell_type
     */
    void removeOneCell();
};

class StochClone: public Clone{
private:
    double drawLogNorm(double mean, double var);
    double drawTruncGamma(double mean, double var);
    double drawTruncDoubleExp(double mean, double var);
protected:
    bool is_mult;
    virtual double setNewBirth(double mean, double var) = 0;
    double drawFromDist(double mean, double var);
    string dist_type;
public:
    virtual void reproduce() = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
    StochClone(CellType& type, bool mult);
    StochClone(CellType& type, double mut, bool mult);
};

class UpdateClone: public Clone{
protected:
    bool is_dead;
    bool has_reproduced;
public:
    UpdateClone(CellType& type);
    bool hasDied(){
        return is_dead;
    }
    bool hasReproduced(){
        return has_reproduced;
    }
    virtual void update(double t) = 0;
    virtual void reproduce() = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
};

class Diffusion1DClone: public UpdateClone{
private:
    double threshold;
    double curr_pos;
    double diffusion;
    double drift;
    bool checkRep(){
        return (diffusion >= 0 && drift >= 0 && curr_pos <= threshold && Clone::checkRep());
    };
public:
    Diffusion1DClone(CellType& type);
    Diffusion1DClone(CellType& type, double b, double mu, double dr, double diff, double thresh, double pos);
    void reproduce();
    void update(double t);
    bool readLine(vector<string>& parsed_line);
};

class EmpiricalClone: public StochClone{
protected:
    double drawEmpirical(double mean, double var);
    bool readDist(string filename);
    virtual double setNewBirth(double mean, double var) = 0;
public:
    virtual void reproduce() = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
    EmpiricalClone(CellType& type, bool mult);
    EmpiricalClone(CellType& type, double mut, bool mult);
};

class SimpleClone: public Clone{
public:
    SimpleClone(CellType& type, double b, double mut, int num_cells);
    SimpleClone(CellType& type);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class TypeSpecificClone: public StochClone{
    // equivalent to a HerResetClone with 1 generation alteration lifetime.
private:
    double mean;
    double var;
protected:
    double setNewBirth(double mean, double var);
public:
    TypeSpecificClone(CellType& type, double mu, double sig, double mut, bool mult);
    TypeSpecificClone(CellType& type, double mu, double sig, double mut, double offset, bool mult);
    TypeSpecificClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HeritableClone: public StochClone{ 
protected:
    double mean;
    double var;
    double setNewBirth(double mean, double var);
public:
    HeritableClone(CellType& type, double mu, double sig, double mut, bool mult, string dist);
    HeritableClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, string dist);
    HeritableClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HerPoissonClone: public HeritableClone{
    // draws a poisson distributed number of alterations per generation
protected:
    double accum_rate;
    double add_alteration();
    
    // adds a Poisson-distributed number of alterations
    double add_alterations();
public:
    HerPoissonClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, double accum, string dist);
    HerPoissonClone(CellType& type, double mu, double sig, double mut, bool mult, double accum, string dist);
    HerPoissonClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HerResetClone: public HeritableClone{
    // draws one new fitness alteration per generation, removes fitness alterations after exactly a given number of generations
private:
    // FIFO queue
    queue<double> active_diff;
    int num_gen_persist;
    
    // called in every reproduction to remove last alteration.
    double reset();
    bool checkRep(){
        return active_diff.size() == num_gen_persist;
    };
public:
    HerResetClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, int num_gen, queue<double>& diffs, string dist);
    HerResetClone(CellType& type, double mu, double sig, double mut, bool mult, int num_gen, queue<double>& diffs, string dist);
    HerResetClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HerResetExpClone: public HerPoissonClone{
    // draws one new fitness alteration per generation, removes fitness alterations after an exponentially-distributed number of generations
private:
    vector<double> active_diff;
    
    // Given as exponential rate constant.
    // The average lifetime is 1/time_constant.
    double time_constant;
    
    // called in every reproduction to choose alteration to remove and remove alterations.
    // DOES NOT ADD ALTERATIONS! add_alteration() does this.
    void reset();
    
    double add_alteration();

    bool checkRep(){
        return (accum_rate > 0 && time_constant > 0);
    };
public:
    HerResetExpClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, double time, double accum, vector<double>& diffs, string dist);
    HerResetExpClone(CellType& type, double mu, double sig, double mut, bool mult, double time, double accum, vector<double>& diffs, string dist);
    HerResetExpClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class TypeEmpiricClone: public EmpiricalClone{
private:
    double mean;
    double var;
protected:
    double setNewBirth(double mean, double var);
public:
    TypeEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult);
    TypeEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult);
    TypeEmpiricClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HerEmpiricClone: public EmpiricalClone{
protected:
    double setNewBirth(double mean, double var);
    double mean;
    double var;
public:
    HerEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult);
    HerEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult);
    HerEmpiricClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

class HerResetEmpiricClone: public HerEmpiricClone{
private:
    // FIFO queue
    queue<double> active_diff;
    int num_gen_persist;
    double reset();
    bool checkRep(){
        return active_diff.size() == num_gen_persist;
    };
public:
    HerResetEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, int num_gen, queue<double>& diffs);
    HerResetEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult, int num_gen, queue<double>& diffs);
    HerResetEmpiricClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

#endif /* clone_h */
