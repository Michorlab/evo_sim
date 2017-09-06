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
#include "MutationHandler.h"
#include <vector>
#include <sstream>
#include <random>
#include <cstdlib>
#include <chrono>
#include <cstdlib>
using namespace std;

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
    death_var = false;
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
    death_var = false;
}

void CList::clearClones(){
    curr_types.clear();
    for (int i = 0; i < max_types; i++){
        curr_types.push_back(NULL);
    }
}

void CList::refreshSim(){
    deleteList();
    clearClones();
    tot_rate = 0;
    time = 0;
    tot_cell_count = 0;
    num_types = 0;
    root = NULL;
    end_node = NULL;
    root_types.clear();
    death_var = false;
}

void CList::insertCellType(CellType& new_type) {
    if (curr_types[new_type.getIndex()]){
        throw "type space conflict";
    }
    if (end_node){
        end_node->setNext(new_type);
    }
    else{
        root = &new_type;
    }
    curr_types[new_type.getIndex()] = &new_type;
    new_type.setPrev(*end_node);
    end_node = &new_type;
    addCells(new_type.getNumCells(), new_type.getBirthRate());
    new_type.setCloneList(*this);
    num_types++;
}

void CList::deleteList()
{
    CellType *to_delete = root;
    CellType *next = to_delete->getNext();
    while(next) {
        to_delete = next;
        next = to_delete->getNext();
        delete to_delete;
    }
    
}

void CList::advance()
{
    uniform_real_distribution<double> runif;
    mut_model->reset();
    double total_death = getTotalDeath();
    if (tot_cell_count == 0){
        tot_rate = 0;
    }
    time += -log(runif(*eng))/(tot_rate + total_death);
    double b_or_d = runif(*eng)*(tot_rate + total_death);
    if (b_or_d < (total_death)){
        if (death_var){
             Clone& dead = chooseDeadVar(total_death);
            if (dead.isSingleCell()){
                delete &dead;
            }
            else{
                dead.removeOneCell();
            }
        }
        else{
             Clone& dead = chooseDead();
            if (dead.isSingleCell()){
                delete &dead;
            }
            else{
                dead.removeOneCell();
            }
        }

    }
    else{
        Clone& mother = chooseReproducer();
        mother.reproduce();
    }
}

Clone& CList::chooseReproducer(){
    uniform_real_distribution<double> runif;
    
    double ran = runif(*eng) * tot_rate;
    CellType *rep_type = root;
    while (rep_type->getNumCells() == 0){
        rep_type = rep_type->getNext();
    }
    Clone *reproducer = rep_type->getRoot();
    double curr_rate = reproducer->getTotalBirth();
    while (curr_rate < ran && reproducer->getNextClone()){
        reproducer = reproducer->getNextClone();
        curr_rate += reproducer->getTotalBirth();
    }
    return *reproducer;
}

double CList::getTotalDeath(){
    if (death_var){
        CellType *curr_type = root;
        double total_d = 0;
        while (curr_type){
            total_d += curr_type->getDeathRate() * curr_type->getNumCells();
            curr_type = curr_type->getNext();
        }
        return total_d;
    }
    else{
        return d*tot_cell_count;
    }
}

Clone& CList::chooseDeadVar(double total_death){
    uniform_real_distribution<double> runif;
    double ran = runif(*eng);
    CellType *dead_type = root;
    while (dead_type->getNumCells() == 0){
        dead_type = dead_type->getNext();
    }
    ran = ran * total_death;
    Clone *dead = dead_type->getRoot();
    double curr_rate = dead->getCellCount() * dead->getDeathRate();
    while (curr_rate < ran && dead->getNextClone()){
        dead = dead->getNextClone();
        curr_rate += dead->getCellCount() * dead->getDeathRate();
    }
    return *dead;

}

Clone& CList::chooseDead(){
    uniform_real_distribution<double> runif;
    double ran = runif(*eng);
    CellType *dead_type = root;
    while (dead_type->getNumCells() == 0){
        dead_type = dead_type->getNext();
    }
    ran = ran * tot_cell_count;
    
    Clone *dead = dead_type->getRoot();
    double curr_rate = dead->getCellCount();
    while (curr_rate < ran && dead->getNextClone()){
        dead = dead->getNextClone();
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
    for (int i=0; i<int(children.size()); i++){
        outfile << children[i]->getIndex() << "\t";
    }
    outfile << endl;
    for (int i=0; i<int(children.size()); i++){
        walkTypesAndWrite(outfile, *children[i]);
    }
}

bool CList::handle_line(vector<string>& parsed_line){
    if (parsed_line[0] == "death"){
        d =stod(parsed_line[1]);
    }
    else if(parsed_line[0] == "max_types"){
        max_types =stoi(parsed_line[1]) + 1;
        clearClones();
    }
    else if(parsed_line[0] == "death_var"){
        death_var = true;
        int type = stoi(parsed_line[1]);
        double death = stod(parsed_line[2]);
        getTypeByIndex(type)->setDeathRate(death);
    }
    else{
        return false;
    }
    
    return true;
}

bool CList::checkInit(){
    return (max_types && mut_model);
}

bool CList::isOneType(){
    return root->getNext();
}

int CList::newestType(){
    return end_node->getIndex();
}

void MoranPop::advance(){
    Clone& dead = chooseDead();
    if (dead.isSingleCell()){
        delete &dead;
    }
    else{
        dead.removeOneCell();
    }
    Clone& mother = chooseReproducer();
    mother.reproduce();
    time++;
}

MoranPop::MoranPop() : CList(){}
