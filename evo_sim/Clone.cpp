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
        cell_type->setRoot(*next_node);
    }
    if (!next_node){
        cell_type->setEnd(*prev_node);
    }
    if (prev_node){
        prev_node->setNext(next_node);
    }
    if (next_node){
        next_node->setPrev(prev_node);
    }
}

Clone* Clone::getNextClone(){
    if (!next_node){
        if (!cell_type->getNext()){
            return NULL;
        }
        else{
            return (cell_type->getNext())->getRoot();
        }
    }
    return next_node;
}

SimpleClone::SimpleClone(CellType& type, double b, double mut, int num_cells) : Clone(type, mut){
    birth_rate = b;
    cell_count = num_cells;
}

Clone::Clone(CellType& type){
    cell_count = 0;
    cell_type = &type;
    mut_prob = 0;
    next_node = NULL;
    prev_node = NULL;
}

SimpleClone::SimpleClone(CellType& type) : Clone(type){};

StochClone::StochClone(CellType& type) : Clone(type){};

TypeSpecificClone::TypeSpecificClone(CellType& type) : StochClone(type){
    mean = 0;
    var = 0;
}

HeritableClone::HeritableClone(CellType& type) : StochClone(type){
    mean = 0;
    var = 0;
}

void SimpleClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(eng3) < mut_prob){
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, birth_rate, mut_prob);
        SimpleClone *new_node = new SimpleClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), mut_handle.getNewMutProb(), 1);
        mut_handle.getNewType().insertClone(*new_node);
    }
    else{
        addCells(1);
    }
}

Clone::Clone(CellType& type, double mut){
    cell_count = 0;
    cell_type = &type;
    mut_prob = mut;
    next_node = NULL;
    prev_node = NULL;
}

StochClone::StochClone(CellType& type, double mut) : Clone(type, mut){}

HeritableClone::HeritableClone(CellType& type, double mu, double sig, double mut) : StochClone(type, mut){
    mean = mu;
    var = sig;
    birth_rate = drawLogNorm(mean, var);
    cell_count = 1;
}

TypeSpecificClone::TypeSpecificClone(CellType& type, double mu, double sig, double mut) : StochClone(type, mut){
    mean = mu;
    var = sig;
    birth_rate = drawLogNorm(mean, var);
    cell_count = 1;
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
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, mean, mut_prob);
        TypeSpecificClone *new_node = new TypeSpecificClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb());
        mut_handle.getNewType().insertClone(*new_node);
    }
    else{
        TypeSpecificClone *new_node = new TypeSpecificClone(*cell_type, mean, var, mut_prob);
        cell_type->insertClone(*new_node);
    }
}

void HeritableClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(eng3) < mut_prob){
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, birth_rate, mut_prob);
        HeritableClone *new_node = new HeritableClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb());
        mut_handle.getNewType().insertClone(*new_node);
    }
    else{
        HeritableClone *new_node = new HeritableClone(*cell_type, birth_rate, var, mut_prob);
        cell_type->insertClone(*new_node);
    }
}

bool SimpleClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone SimpleClone [type_id] [num_cells] [birth_rate] [mut_rate]
    try{
        cell_count = stoi(parsed_line[0]);
        birth_rate = stod(parsed_line[1]);
        mut_prob = stod(parsed_line[2]);
    }
    catch (const std::invalid_argument&){
        return false;
    }
    catch (const std::out_of_range&){
        return false;
    }
    return checkRep();
}

bool TypeSpecificClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone TypeSpecificClone [type_id] [num_cells] [mean] [var] [mut_rate]
    cell_count = 1;
    try{
        mean = stod(parsed_line[1]);
        var = stod(parsed_line[2]);
        mut_prob = stod(parsed_line[3]);
    }
    catch (const std::invalid_argument&){
        return false;
    }
    catch (const std::out_of_range&){
        return false;
    }
    birth_rate = drawLogNorm(mean, var);
    return checkRep();
}

bool HeritableClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HeritableClone [type_id] [num_cells] [mean] [var] [mut_rate]
    cell_count = 1;
    try{
        mean = stod(parsed_line[1]);
        var = stod(parsed_line[2]);
        mut_prob = stod(parsed_line[3]);
    }
    catch (const std::invalid_argument&){
        return false;
    }
    catch (const std::out_of_range&){
        return false;
    }
    birth_rate = drawLogNorm(mean, var);
    return checkRep();
}
