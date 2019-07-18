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
    recalc_birth = false;
    prev_fit = 0;
    new_fit = 0;
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
    recalc_birth = false;
    prev_fit = 0;
    new_fit = 0;
    new_type = 0;
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
    prev_fit = 0;
    new_fit = 0;
    new_type = 0;
}

void SexReprPop::refreshSim(){
    CList::refreshSim();
    is_extinct = false;
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

double CList::getTotalBirth(){
    if (recalc_birth){
        CellType *rep_type = root;
        while (rep_type->getNumCells() == 0){
            rep_type = rep_type->getNext();
        }
        Clone *reproducer = rep_type->getRoot();
        double curr_rate = reproducer->getTotalBirth();
        while (reproducer->getNextClone()){
            reproducer = reproducer->getNextClone();
            curr_rate += reproducer->getTotalBirth();
        }
        return curr_rate;
    }
    else{
        return tot_rate;
    }
}

void CList::killCell(Clone& dead){
    if (dead.isSingleCell()){
        delete &dead;
    }
    else{
        dead.removeOneCell();
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
    double tot_birth = getTotalBirth();
    time += -log(runif(*eng))/(tot_birth + total_death);
    double b_or_d = runif(*eng)*(tot_birth + total_death);
    if (b_or_d < (total_death)){
        if (death_var){
            Clone& dead = chooseDeadVar(total_death);
            killCell(dead);

        }
        else{
            Clone& dead = chooseDead();
            killCell(dead);
        }

    }
    else{
        Clone& mother = chooseReproducer();
        prev_fit = mother.getBirthRate();
        mother.reproduce();
        new_fit = mother.getBirthRate();
        if (mut_model->has_mut()){
            new_type = mut_model->getNewType().getIndex();
        }
    }
}

Clone& CList::chooseReproducer(){
    uniform_real_distribution<double> runif;
    
    double ran = runif(*eng) * getTotalBirth();
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

/*
void CList::walkTypesAndWrite(ofstream& outfile, CellType& root){
    outfile << root.getIndex() << ", " << root.isExtinct() << ", ";
    std::vector<CellType *> children = root.getChildren();
    for (int i=0; i<int(children.size()); i++){
        if (children[i]->getParent()->getIndex() == root.getIndex()){
            outfile << children[i]->getIndex() << ", ";
        }
    }
    outfile << endl;
    for (int i=0; i<int(children.size()); i++){
        if (children[i]->getParent()->getIndex() == root.getIndex()){
            walkTypesAndWrite(outfile, *children[i]);
        }
    }
}
*/

void CList::walkTypesAndWrite(ofstream& outfile, CellType& root){
    for (int i=0; i<max_types; i++){
        if (hasCellType(i)){
            outfile << i << ", " << getTypeByIndex(i)->getNumCells() << ", ";
            if (getTypeByIndex(i)->getParent()){
                outfile << getTypeByIndex(i)->getParent()->getIndex() << endl;
            }
            else{
                outfile << endl;
            }
        }
    }
}
 
bool CList::handle_line(vector<string>& parsed_line){
    if (parsed_line[0] == "death"){
        d =stod(parsed_line[1]);
    }
    else if (parsed_line[0] == "recalc_birth"){
        recalc_birth = true;
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

bool UpdateAllPop::checkInit(){
    return (CList::checkInit() && timestep_length);
}

bool CList::isOneType(){
    return root->getNext();
}

int CList::newestType(){
    return end_node->getIndex();
}

void MoranPop::advance(){
    mut_model->reset();
    Clone& dead = chooseDead();
    killCell(dead);
    Clone& mother = chooseReproducer();
    prev_fit = mother.getBirthRate();
    mother.reproduce();
    new_fit = mother.getBirthRate();
    if (mut_model->has_mut()){
        new_type = mut_model->getNewType().getIndex();
    }
    time++;
}

MoranPop::MoranPop() : CList(){}

SexReprPop::SexReprPop() : CList(){
    std::vector<int> male_types = std::vector<int>();
    std::vector<int> female_types = std::vector<int>();
    is_extinct = false;
}

UpdateAllPop::UpdateAllPop() : CList(){
    timestep_length = 0;
}

void UpdateAllPop::advance(){
    mut_model->reset();
    std::vector<Clone *> reproducers = std::vector<Clone *>();
    std::vector<Clone *> dead = std::vector<Clone *>();
    
    CellType *curr_type = root;
    while (curr_type->getNumCells() == 0){
        curr_type = curr_type->getNext();
    }
    Clone *curr = curr_type->getRoot();
    while (curr->getNextClone()){
        curr->update(timestep_length);
        if (curr->hasDied()){
            dead.push_back(curr);
        }
        else if (curr->hasReproduced()){
            reproducers.push_back(curr);
        }
        curr = curr->getNextClone();
    }
    curr->update(timestep_length);
    if (curr->hasDied()){
        dead.push_back(curr);
    }
    else if (curr->hasReproduced()){
        reproducers.push_back(curr);
    }
    
    for (int i=0; i<reproducers.size(); i++){
        reproducers.back()->reproduce();
        reproducers.pop_back();
    }
    
    for (int i=0; i<dead.size(); i++){
        killCell(*dead.back());
        dead.pop_back();
    }
    
    time += timestep_length;
}

bool UpdateAllPop::handle_line(vector<string>& parsed_line){
    if (parsed_line[0] == "timestep"){
        timestep_length =stod(parsed_line[1]);
    }
    else{
        return CList::handle_line(parsed_line);
    }
    
    return true;
}

bool SexReprPop::handle_line(vector<string>& parsed_line){
    if (parsed_line[0] == "male_types"){
        for (int i=1; i<parsed_line.size(); i++){
            male_types.push_back(stoi(parsed_line[i]));
        }
    }
    else if (parsed_line[0] == "female_types"){
        for (int i=1; i<parsed_line.size(); i++){
            female_types.push_back(stoi(parsed_line[i]));
        }
    }
    else{
        return CList::handle_line(parsed_line);
    }
    
    return true;
}

void SexReprPop::advance(){
    mut_model->reset();
    std::vector<SexReprClone *> new_cells = std::vector<SexReprClone *>();
    std::vector<int> type_indices = std::vector<int>();
    for (int i=0; i<tot_cell_count; i++){
        SexReprClone* mother = &chooseMother();
        SexReprClone* father = &chooseFather();
        SexReprClone& new_cell = mother->reproduce(*father);
        new_cells.push_back(&new_cell);
        type_indices.push_back(new_cell.getType().getIndex());
    }
    double prev_time = time;
    refreshSim();
    for (vector<int>::iterator it = male_types.begin(); it != male_types.end(); ++it){
        CellType *new_type = new CellType(*it, NULL);
        insertCellType(*new_type);
    }
    for (vector<int>::iterator it = female_types.begin(); it != female_types.end(); ++it){
        CellType *new_type = new CellType(*it, NULL);
        insertCellType(*new_type);
    }
    for (int i=0; i<new_cells.size(); i++){
        SexReprClone* new_cell = new_cells[i];
        int index = type_indices[i];
        CellType *new_type = getTypeByIndex(index);
        new_cell->setType(*new_type);
        new_cell->getType().insertClone(*new_cell);
    }
    bool males_extinct = true;
    bool females_extinct = true;
    for (vector<int>::iterator it = male_types.begin(); it != male_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        males_extinct = males_extinct && curr_type->isExtinct();
    }
    for (vector<int>::iterator it = female_types.begin(); it != female_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        females_extinct = females_extinct && curr_type->isExtinct();
    }
    is_extinct = males_extinct && females_extinct;
    time = prev_time + 1;
    recalc_birth = true;
}

bool SexReprPop::checkInit(){
    bool males_extinct = false;
    bool females_extinct = false;
    for (vector<int>::iterator it = male_types.begin(); it != male_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        males_extinct = males_extinct || curr_type->isExtinct();
    }
    for (vector<int>::iterator it = female_types.begin(); it != female_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        females_extinct = females_extinct || curr_type->isExtinct();
    }
    is_extinct = males_extinct || females_extinct;
    return !is_extinct && CList::checkInit();
}

SexReprClone& SexReprPop::chooseReproducerVector(vector<int> possible_types){
    uniform_real_distribution<double> runif;
    double total_birth_vect = 0;
    for (vector<int>::iterator it = possible_types.begin(); it != possible_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        if (!curr_type || curr_type->isExtinct()){
            continue;
        }
        total_birth_vect += curr_type->getBirthRate();
    }
    double ran = runif(*eng) * total_birth_vect;
    double curr_rate = 0;
    SexReprClone *reproducer;
    for (vector<int>::iterator it = possible_types.begin(); it != possible_types.end(); ++it){
        CellType* curr_type = getTypeByIndex(*it);
        if (!curr_type || curr_type->isExtinct()){
            continue;
        }
        reproducer = (SexReprClone*)curr_type->getRoot();
        curr_rate += reproducer->getTotalBirth();
        while (curr_rate < ran && reproducer->getNextClone()){
            if (reproducer->getNextClone()->getType().getIndex() != reproducer->getType().getIndex()){
                break;
            }
            reproducer = (SexReprClone*)reproducer->getNextClone();
            curr_rate += reproducer->getTotalBirth();
        }
        if (curr_rate > ran){
            break;
        }
    }
    return *reproducer;
}

SexReprClone& SexReprPop::chooseFather(){
    return chooseReproducerVector(male_types);
}

SexReprClone& SexReprPop::chooseMother(){
    return chooseReproducerVector(female_types);
}

void SexReprPop::addMaleType(int type_index){
    if(male_types.size() == 0 || std::find(male_types.begin(), male_types.end(), type_index) == male_types.end()) {
        male_types.push_back(type_index);
    }
    return;
}

void SexReprPop::addFemaleType(int type_index){
    if(female_types.size() == 0 || std::find(female_types.begin(), female_types.end(), type_index) == female_types.end()) {
        female_types.push_back(type_index);
    }
    return;
}
