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

using namespace std;

class CList;
class CellType;
class MutationHandler;

class Clone{
protected:
    long long cell_count;
    CellType *cell_type;
    
    double birth_rate;
    double mut_prob;
    
    CList *clone_list;
    Clone *next_node;
    Clone *prev_node;
    
public:
    /* adjusts the linked list so the clone is no longer in the list. adjusts root and end of clone_list as appropriate.
     MODIFIES clone_list.
     */
    ~Clone();
    
    Clone(CellType& type, double mut, CList& pop);
    
    virtual void reproduce() = 0;
    
    /* removes one cell from this clone's population
     should not be called if there is <=1 cell left in the clone
     MODIFIES clone_list
     */
    void removeOneCell();
    
    double getBirthRate(){
        return birth_rate;
    }
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
    Clone& getNext(){
        return *next_node;
    }
    Clone& getPrev(){
        return *prev_node;
    }
    CellType& getType(){
        return *cell_type;
    }
    void addCells(int num_cells);
};

class StochClone: public Clone{
protected:
    double drawLogNorm(double mean, double var);
public:
    virtual void reproduce() = 0;
    StochClone(CellType& type, double mut, CList& pop);
};

class SimpleClone: public Clone{
public:
    SimpleClone(CellType& type, double b, double mut, int num_cells, CList& pop);
    void reproduce();
};

class TypeSpecificClone: public StochClone{
private:
    double mean;
    double var;
public:
    TypeSpecificClone(CellType& type, double mu, double sig, double mut, CList& pop);
    void reproduce();
};

class HeritableClone: public StochClone{
private:
    double mean;
    double var;
public:
    HeritableClone(CellType& type, double mu, double sig, double mut, CList& pop);
    void reproduce();
};

#endif /* clone_h */
