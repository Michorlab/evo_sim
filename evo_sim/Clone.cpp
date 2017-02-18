//
//  Clone.cpp
//  variable fitness branching process
//
//  Created by Debra Van Egeren on 2/17/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#include "Clone.h"
#include "CList.h"
#include "main.h"
#include <random>
#include <vector>
#include <chrono>
#include <cstdlib>
using namespace std;

int seed3 =  std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 eng3(seed3);

void Clone::removeOneCell(){
    if (cell_count <= 1){
        throw "death error";
    }
    else{
        cell_count--;
        cell_type->subtractOneCell(birth_rate);
    }
}

Clone::~Clone(){
    cell_type->subtractOneCell(birth_rate);
    if (!prev_node){
        clone_list->setRoot(next_node);
    }
    if (!next_node){
        clone_list->setEnd(prev_node);
    }
    if (prev_node){
        prev_node->setNext(next_node);
    }
    if (next_node){
        next_node->setPrev(prev_node);
    }
}

SimpleClone::SimpleClone(CellType& type, double b, double mut, int num_cells, CList& pop) : Clone(type, mut, pop){
    cell_type = &type;
    birth_rate = b;
    next_node = NULL;
    prev_node = &pop.getEnd();
}

void SimpleClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(eng3) < mut_prob){
        MutationHandler& mut_handle = clone_list->getMutHandler();
        mut_handle.generateMutant(*cell_type, *clone_list, birth_rate, mut_prob);
        SimpleClone *new_node = new SimpleClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), mut_handle.getNewMutProb(), 1, *clone_list);
        clone_list->insertNode(*new_node);
    }
    else{
        addCells(1);
    }
}

Clone::Clone(CellType& type, double mut, CList& pop){
    cell_count = 0;
    cell_type = &type;
    mut_prob = mut;
    clone_list = &pop;
    next_node = NULL;
    prev_node = &pop.getEnd();
}

StochClone::StochClone(CellType& type, double mut, CList& pop) : Clone(type, mut, pop){}

HeritableClone::HeritableClone(CellType& type, double mu, double sig, double mut, CList& pop) : StochClone(type, mut, pop){
    mean = mu;
    var = sig;
    birth_rate = drawLogNorm(mean, var);
    addCells(1);
}

TypeSpecificClone::TypeSpecificClone(CellType& type, double mu, double sig, double mut, CList& pop) : StochClone(type, mut, pop){
    mean = mu;
    var = sig;
    birth_rate = drawLogNorm(mean, var);
    addCells(1);
}

void Clone::addCells(int num_cells){
    cell_count+=num_cells;
    cell_type->addCells(birth_rate, num_cells);
}

double StochClone::drawLogNorm(double mean, double var){
    double loc = log(pow(mean, 2.0)/sqrt(var+pow(mean,2.0)));
    double scale = log(1.0+var/pow(mean,2.0));
    normal_distribution<double> norm(loc,scale);
    return exp(norm(eng3));
}

void TypeSpecificClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(eng3) < mut_prob){
        MutationHandler& mut_handle = clone_list->getMutHandler();
        mut_handle.generateMutant(*cell_type, *clone_list, mean, mut_prob);
        TypeSpecificClone *new_node = new TypeSpecificClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), *clone_list);
        clone_list->insertNode(*new_node);
    }
    else{
        TypeSpecificClone *new_node = new TypeSpecificClone(*cell_type, mean, var, mut_prob, *clone_list);
        clone_list->insertNode(*new_node);
    }
}

void HeritableClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(eng3) < mut_prob){
        MutationHandler& mut_handle = clone_list->getMutHandler();
        mut_handle.generateMutant(*cell_type, *clone_list, birth_rate, mut_prob);
        HeritableClone *new_node = new HeritableClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), *clone_list);
        clone_list->insertNode(*new_node);
    }
    else{
        TypeSpecificClone *new_node = new TypeSpecificClone(*cell_type, birth_rate, var, mut_prob, *clone_list);
        clone_list->insertNode(*new_node);
    }
}
