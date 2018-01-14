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
protected:
    bool is_mult;
    double drawLogNorm(double mean, double var);
    virtual double setNewBirth(double mean, double var) = 0;
public:
    virtual void reproduce() = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
    StochClone(CellType& type, bool mult);
    StochClone(CellType& type, double mut, bool mult);
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
private:
    double mean;
    double var;
protected:
    double setNewBirth(double mean, double var);
public:
    HeritableClone(CellType& type, double mu, double sig, double mut, bool mult);
    HeritableClone(CellType& type, double mu, double sig, double mut, double offset, bool mult);
    HeritableClone(CellType& type, bool mult);
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
private:
    double mean;
    double var;
protected:
    double setNewBirth(double mean, double var);
public:
    HerEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult);
    HerEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult);
    HerEmpiricClone(CellType& type, bool mult);
    void reproduce();
    bool readLine(vector<string>& parsed_line);
};

#endif /* clone_h */
