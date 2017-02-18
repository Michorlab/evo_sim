//
//  CList.cpp
//  variable fitness branching process
//
//  Created by Debra Van Egeren on 2/17/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#include "CList.h"
#include "Clone.h"
#include "main.h"
#include <vector>
#include <sstream>
#include <random>
#include <cstdlib>
#include <chrono>
#include <cstdlib>
using namespace std;

int seed2 =  std::chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 eng2(seed2);

CList::CList(double death, MutationHandler& mut_handle, int max){
    d = death;
    tot_rate = 0;
    time = 0;
    max_types = max;
    num_types = 0;
    curr_types.reserve(max);
    clearClones();
    tot_cell_count = 0;
    mut_model = &mut_handle;
    root = NULL;
    end_node = NULL;
}

CList::CList(){
    tot_rate = 0;
    time = 0;
    num_types = 0;
    tot_cell_count = 0;
    root = NULL;
    end_node = NULL;
    mut_model = NULL;
    d = 0;
}

void CList::clearClones(){
    for (int i = 0; i < max_types; i++){
        delete curr_types[i];
        curr_types[i] = NULL;
    }
}

void CList::refreshSim(){
    deleteList();
    clearClones();
    tot_rate = 0;
    time = 0;
    tot_cell_count = 0;
    num_types = 0;
}

void CList::insertCellType(CellType& new_type) {
    if (curr_types[new_type.getIndex()]){
        throw "type space conflict";
    }
    end_node->setNext(new_type);
    new_type.setPrev(*end_node);
    end_node = &new_type;
    addCells(new_type.getNumCells(), new_type.getBirthRate());
}

void CList::deleteList()
{
    CellType *pnode = root;
    while(root!=NULL) {
        pnode = root;
        root = &root->getNext();
        delete pnode;
    }
    
}

void CList::advance()
{
    uniform_real_distribution<double> runif;
    
    time += -log(runif(eng2))/tot_rate;
    double b_or_d = runif(eng2)*tot_rate;
    if (b_or_d < d * tot_cell_count){
        Clone& dead = chooseDead();
        if (dead.isSingleCell()){
            delete &dead;
        }
        else{
            dead.removeOneCell();
        }
    }
    else{
        Clone& mother = chooseReproducer();
        mother.reproduce();
    }
}

Clone& CList::chooseReproducer(){
    uniform_real_distribution<double> runif;
    
    double ran = runif(eng2) * tot_rate;
    Clone *reproducer = &root->getRoot();
    double curr_rate = reproducer->getTotalBirth();
    while (curr_rate < ran && reproducer){
        reproducer = &reproducer->getNextClone();
        curr_rate += reproducer->getTotalBirth();
    }
    return *reproducer;
}

Clone& CList::chooseDead(){
    uniform_real_distribution<double> runif;
    
    double ran = runif(eng2) * tot_cell_count;
    Clone *dead = &root->getRoot();
    double curr_rate = dead->getCellCount();
    while (curr_rate < ran && dead){
        dead = &dead->getNextClone();
        curr_rate += dead->getCellCount();
    }
    return *dead;
    
}

int CList::getNextType(){
    for (int i=0; i < max_types; i++){
        if (!curr_types[i]){
            return i;
        }
    }
    throw "tried to get a new type at max_types";
}

void CList::addCells(int num_cells, double b){
    tot_rate += b * num_cells;
    tot_cell_count += num_cells;
}

void CList::removeCell(double b){
    tot_rate -= b;
    tot_cell_count --;
}

void CList::walkTypesAndWrite(ofstream& outfile, CellType& root){
    outfile << root.getIndex() << "\t" << root.isExtinct() << "\t";
    std::vector<CellType *> children = root.getChildren();
    for (int i=0; i<children.size(); i++){
        outfile << children[i]->getIndex() << "\t";
    }
    outfile << endl;
    for (int i=0; i<children.size(); i++){
        walkTypesAndWrite(outfile, *children[i]);
    }
}

bool CList::handle_line(vector<string>& parsed_line){
    if (parsed_line[0] == "death"){
        d = stod(parsed_line[1]);
    }
    else if(parsed_line[0] == "max_types"){
        max_types = stoi(parsed_line[1]) + 1;
        clearClones();
    }
    else{
        return false;
    }
    return true;
}

bool CList::checkInit(){
    return (max_types && mut_model);
}
