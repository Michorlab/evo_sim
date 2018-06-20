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
#include <queue>
#include <chrono>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include "MutationHandler.h"
using namespace std;

void Clone::removeOneCell(){
    cell_count--;
    cell_type->subtractOneCell(birth_rate);
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
            CellType* next_type = cell_type->getNext();
            Clone* next_clone = next_type->getRoot();
            while (next_type && !next_clone){
                next_type = next_type->getNext();
                next_clone = next_type->getRoot();
            }
            return next_clone;
        }
    }
    return next_node;
}

double Clone::getDeathRate(){
    return cell_type->getDeathRate();
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

StochClone::StochClone(CellType& type, bool mult) : Clone(type){
    is_mult = mult;
    dist_type = "lognorm";
};

EmpiricalClone::EmpiricalClone(CellType& type, bool mult) : StochClone(type, mult){};

TypeSpecificClone::TypeSpecificClone(CellType& type, bool mult) : StochClone(type, mult){
    mean = 0;
    var = 0;
}

HeritableClone::HeritableClone(CellType& type, bool mult) : StochClone(type, mult){
    mean = 0;
    var = 0;
}

HerResetClone::HerResetClone(CellType& type, bool mult) : HeritableClone(type, mult){
    num_gen_persist = 0;
    active_diff = queue<double>();
}

HerResetExpClone::HerResetExpClone(CellType& type, bool mult) : HeritableClone(type, mult){
    time_constant = 0;
    accum_rate = 0;
    active_diff = vector<double>();
}

HerResetEmpiricClone::HerResetEmpiricClone(CellType& type, bool mult) : HerEmpiricClone(type, mult){
    num_gen_persist = 0;
    active_diff = queue<double>();
}

TypeEmpiricClone::TypeEmpiricClone(CellType& type, bool mult) : EmpiricalClone(type, mult){
    mean = 0;
    var = 0;
}

HerEmpiricClone::HerEmpiricClone(CellType& type, bool mult) : EmpiricalClone(type, mult){
    mean = 0;
    var = 0;
}

void SimpleClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, birth_rate, mut_prob);
        if (mut_handle.getNewType().getEnd() && mut_handle.getNewType().getEnd()->getBirthRate()==mut_handle.getNewBirthRate() && mut_handle.getNewType().getEnd()->getMutProb()== mut_handle.getNewMutProb()){
            mut_handle.getNewType().getEnd()->addCells(1);
        }
        else{
            SimpleClone *new_node = new SimpleClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), mut_handle.getNewMutProb(), 1);
            mut_handle.getNewType().insertClone(*new_node);
        }
        
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

StochClone::StochClone(CellType& type, double mut, bool mult) : Clone(type, mut){
    is_mult = mult;
    dist_type = "lognorm";
}

EmpiricalClone::EmpiricalClone(CellType& type, double mut, bool mult) : StochClone(type, mut, mult){};

HeritableClone::HeritableClone(CellType& type, double mu, double sig, double mut, bool mult, string dist) : StochClone(type, mut, mult){
    dist_type = dist;
    mean = mu;
    var = sig;
    setNewBirth(mean, var);
    cell_count = 1;
}

HerEmpiricClone::HerEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult) : EmpiricalClone(type, mut, mult){
    mean = mu;
    var = sig;
    setNewBirth(mean, var);
    cell_count = 1;
}

HeritableClone::HeritableClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, string dist) : StochClone(type, mut, mult){
    mean = mu;
    var = sig;
    if (is_mult){
        birth_rate = mean * offset;
    }
    else{
        birth_rate = mean + offset;
    }
    cell_count = 1;
    dist_type = dist;
}

HerResetClone::HerResetClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, int num_gen, queue<double>& diffs, string dist) : HeritableClone(type, mu, sig, mut, offset, mult, dist){
    num_gen_persist = num_gen;
    active_diff = queue<double>(diffs);
    if (!HerResetClone::checkRep()){
        throw "mismanaged reset queue";
    }
}

HerResetExpClone::HerResetExpClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, double time, double accum, vector<double>& diffs, string dist) : HeritableClone(type, mu, sig, mut, offset, mult, dist){
    time_constant = time;
    accum_rate = accum;
    active_diff = vector<double>(diffs);
    if (!HerResetExpClone::checkRep()){
        throw "bad time constant for HerResetExp";
    }
}

HerResetEmpiricClone::HerResetEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult, int num_gen, queue<double>& diffs) : HerEmpiricClone(type, mu, sig, mut, offset, mult){
    num_gen_persist = num_gen;
    active_diff = queue<double>(diffs);
    if (!HerResetEmpiricClone::checkRep()){
        throw "mismanaged reset queue";
    }
}

HerResetClone::HerResetClone(CellType& type, double mu, double sig, double mut, bool mult, int num_gen, queue<double>& diffs, string dist) : HeritableClone(type, mult){
    mean = mu;
    var = sig;
    mut_prob = mut;
    birth_rate = mu;
    num_gen_persist = num_gen;
    active_diff = queue<double>(diffs);
    dist_type = dist;
    cell_count = 1;
    if (!HerResetClone::checkRep()){
        throw "mismanaged reset queue";
    }
}

HerResetExpClone::HerResetExpClone(CellType& type, double mu, double sig, double mut, bool mult, double time, double accum, vector<double>& diffs, string dist) : HeritableClone(type, mult){
    mean = mu;
    var = sig;
    mut_prob = mut;
    birth_rate = mu;
    accum_rate = accum;
    time_constant = time;
    active_diff = vector<double>(diffs);
    dist_type = dist;
    cell_count = 1;
    if (!HerResetExpClone::checkRep()){
        throw "bad time constant for HerResetExp";
    }
}

HerResetEmpiricClone::HerResetEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult, int num_gen, queue<double>& diffs) : HerEmpiricClone(type, mult){
    mean = mu;
    var = sig;
    mut_prob = mut;
    birth_rate = mu;
    num_gen_persist = num_gen;
    active_diff = queue<double>(diffs);
    cell_count = 1;
    if (!HerResetEmpiricClone::checkRep()){
        throw "mismanaged reset queue";
    }
}

HerEmpiricClone::HerEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult) : EmpiricalClone(type, mut, mult){
    mean = mu;
    var = sig;
    birth_rate = mean * offset;
    cell_count = 1;
}

TypeSpecificClone::TypeSpecificClone(CellType& type, double mu, double sig, double mut, bool mult) : StochClone(type, mut, mult){
    mean = mu;
    var = sig;
    setNewBirth(mean, var);
    cell_count = 1;
}

TypeEmpiricClone::TypeEmpiricClone(CellType& type, double mu, double sig, double mut, bool mult) : EmpiricalClone(type, mut, mult){
    mean = mu;
    var = sig;
    setNewBirth(mean, var);
    cell_count = 1;
}

TypeSpecificClone::TypeSpecificClone(CellType& type, double mu, double sig, double mut, double offset, bool mult) : StochClone(type, mut, mult){
    mean = mu;
    var = sig;
    birth_rate = mean * offset;
    cell_count = 1;
}

TypeEmpiricClone::TypeEmpiricClone(CellType& type, double mu, double sig, double mut, double offset, bool mult) : EmpiricalClone(type, mut, mult){
    mean = mu;
    var = sig;
    birth_rate = mean * offset;
    cell_count = 1;
}

void Clone::addCells(int num_cells){
    cell_count+=num_cells;
    cell_type->addCells(num_cells, birth_rate);
}

double StochClone::drawLogNorm(double mean, double var){
    double loc = log(pow(mean, 2.0)/sqrt(var+pow(mean,2.0)));
    double scale = sqrt(log(1.0+var/pow(mean,2.0)));
    normal_distribution<double> norm(loc, scale);
    double to_return = exp(norm(*eng));
    if (to_return < 0){
        return 0;
    }
    return to_return;
}

double StochClone::drawTruncDoubleExp(double mean, double var){
    double lambda = 1.0/sqrt(var/2.0);
    exponential_distribution<double> expo(lambda);
    std::bernoulli_distribution bern(0.5);
    double to_return = expo(*eng);
    if (bern(*eng)){
        to_return = -to_return;
    }
    to_return += mean;
    if (to_return < 0){
        return 0;
    }
    return to_return;
}

double StochClone::drawTruncGamma(double mean, double var){
    double beta = var/mean;
    double alpha = mean/beta;
    std::gamma_distribution<double> gam(alpha,beta);
    double to_return = gam(*eng);
    if (to_return < 0){
        return 0;
    }
    return to_return;
}

void TypeSpecificClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        removeOneCell();
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, mean, mut_prob);
        double offset = setNewBirth(mean, var);
        TypeSpecificClone *new_node = new TypeSpecificClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        TypeSpecificClone *new_node = new TypeSpecificClone(*cell_type, mean, var, mut_prob, is_mult);
        removeOneCell();
        birth_rate = new_node->getBirthRate();
        cell_type->insertClone(*new_node);
        addCells(1);
    }
}

double TypeSpecificClone::setNewBirth(double mean, double var){
    double offset = 0;
    if (is_mult){
        offset = drawFromDist(1, var);
        if (offset < 0){
            offset = 0;
        }
        birth_rate = offset * mean;
    }
    else{
        offset = drawFromDist(0, var);
        if (birth_rate + offset < 0){
            offset = -birth_rate;
        }
        birth_rate = offset + mean;
    }
    return offset;
}

double TypeEmpiricClone::setNewBirth(double mean, double var){
    double offset = 0;
    if (is_mult){
        offset = drawEmpirical(1, var);
        if (offset < 0){
            offset = 0;
        }
        birth_rate = offset * mean;
    }
    else{
        offset = drawEmpirical(0, var);
        if (birth_rate + offset < 0){
            offset = -birth_rate;
        }
        birth_rate = offset + mean;
    }
    return offset;
}

void TypeEmpiricClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        removeOneCell();
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, mean, mut_prob);
        double offset = setNewBirth(mean, var);
        TypeEmpiricClone *new_node = new TypeEmpiricClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        TypeEmpiricClone *new_node = new TypeEmpiricClone(*cell_type, mean, var, mut_prob, is_mult);
        removeOneCell();
        birth_rate = new_node->getBirthRate();
        cell_type->insertClone(*new_node);
        addCells(1);
    }
}

void HeritableClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, birth_rate, mut_prob);
        removeOneCell();
        double offset = setNewBirth(birth_rate, var);
        HeritableClone *new_node = new HeritableClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult, dist_type);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        removeOneCell();
        HeritableClone *new_node = new HeritableClone(*cell_type, birth_rate, var, mut_prob, is_mult, dist_type);
        birth_rate = new_node->getBirthRate();
        addCells(1);
        cell_type->insertClone(*new_node);
    }
}

double HerResetClone::reset(){
    double to_remove = active_diff.front();
    removeOneCell();
    if (is_mult){
        birth_rate = birth_rate/to_remove;
    }
    else{
        birth_rate = birth_rate - to_remove;
    }
    active_diff.pop();
    double offset = setNewBirth(birth_rate, var);
    active_diff.push(offset);
    return offset;
}

void HerResetExpClone::reset(){

    
    double to_remove;
    vector<int> index_to_remove = vector<int>();
    uniform_real_distribution<double> runif;
    
    if (is_mult){
        to_remove = 1;
        for (int i=0; i<active_diff.size(); i++){
            if (runif(*eng) < time_constant){
                to_remove *= active_diff.at(i);
                index_to_remove.push_back(i);
            }
        }
    }
    else{
        to_remove = 0;
        for (int i=0; i<active_diff.size(); i++){
            if (runif(*eng) < time_constant){
                to_remove += active_diff.at(i);
                index_to_remove.push_back(i);
            }
        }
    }
    
    removeOneCell();
    if (index_to_remove.size()==0){
        return;
    }

    
    if (is_mult){
        birth_rate = birth_rate/to_remove;
    }
    else{
        birth_rate = birth_rate - to_remove;
    }
    
    for (int i=index_to_remove.size()-1; i>=0;i--){
        active_diff.erase(active_diff.begin() + i);
    }
}

double HerResetExpClone::add_alteration(){
    double offset = setNewBirth(birth_rate, var);
    active_diff.push_back(offset);
    return offset;
}

double HerResetExpClone::add_alterations(){
    std::poisson_distribution<int> rpois(accum_rate);
    
    double offset;
    if (is_mult){
        offset = 1;
        for (int i=0; i<rpois(*eng); i++){
            offset *= add_alteration();
        }
    }
    else{
        offset=0;
        for (int i=0; i<rpois(*eng); i++){
            offset += add_alteration();
        }
    }
    return offset;
}

double HerResetEmpiricClone::reset(){
    double to_remove = active_diff.front();
    removeOneCell();
    if (is_mult){
        birth_rate = birth_rate/to_remove;
    }
    else{
        birth_rate = birth_rate - to_remove;
    }
    active_diff.pop();
    double offset = setNewBirth(birth_rate, var);
    active_diff.push(offset);
    return offset;
}

void HerResetClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        double offset = reset();
        MutationHandler& mut_handle = cell_type->getMutHandler();
        if (is_mult){
            mut_handle.generateMutant(*cell_type, birth_rate/offset, mut_prob);
        }
        else{
            mut_handle.generateMutant(*cell_type, birth_rate - offset, mut_prob);
        }
        HerResetClone *new_node = new HerResetClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult, num_gen_persist, active_diff, dist_type);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        reset();
        HerResetClone *new_node = new HerResetClone(*cell_type, birth_rate, var, mut_prob, is_mult, num_gen_persist, active_diff, dist_type);
        addCells(1);
        cell_type->insertClone(*new_node);
    }
}

void HerResetExpClone::reproduce(){
    uniform_real_distribution<double> runif;
    
    if (runif(*eng) < mut_prob){
        reset();
        double offset = add_alterations();
        MutationHandler& mut_handle = cell_type->getMutHandler();
        if (is_mult){
            mut_handle.generateMutant(*cell_type, birth_rate/offset, mut_prob);
        }
        else{
            mut_handle.generateMutant(*cell_type, birth_rate - offset, mut_prob);
        }
        HerResetExpClone *new_node = new HerResetExpClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult, time_constant, accum_rate, active_diff, dist_type);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        reset();
        add_alterations();
        HerResetExpClone *new_node = new HerResetExpClone(*cell_type, birth_rate, var, mut_prob, is_mult, time_constant, accum_rate, active_diff, dist_type);
        addCells(1);
        cell_type->insertClone(*new_node);
    }
}

void HerResetEmpiricClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        double offset = reset();
        MutationHandler& mut_handle = cell_type->getMutHandler();
        if (is_mult){
            mut_handle.generateMutant(*cell_type, birth_rate/offset, mut_prob);
        }
        else{
            mut_handle.generateMutant(*cell_type, birth_rate - offset, mut_prob);
        }
        HerResetEmpiricClone *new_node = new HerResetEmpiricClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult, num_gen_persist, active_diff);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        reset();
        HerResetEmpiricClone *new_node = new HerResetEmpiricClone(*cell_type, birth_rate, var, mut_prob, is_mult, num_gen_persist, active_diff);
        addCells(1);
        cell_type->insertClone(*new_node);
    }
}

double StochClone::drawFromDist(double mean, double var){
    if (dist_type=="lognorm"){
        return drawLogNorm(mean, var);
    }
    else if (dist_type=="gamma"){
        return drawTruncGamma(mean, var);
    }
    else if (dist_type=="expo"){
        return drawTruncDoubleExp(mean, var);
    }
    else{
        std::cout << dist_type << std::endl; 
        throw "bad dist type";
        return 0;
    }
}

double HeritableClone::setNewBirth(double mean, double var){
    double offset = 0;
    if (is_mult){
        offset = drawFromDist(1, var);
        if (offset < 0){
            offset = 0;
        }
        birth_rate = offset * mean;
    }
    else{
        offset = drawFromDist(0, var);
        if (birth_rate + offset < 0){
            offset = -birth_rate;
        }
        birth_rate = offset + mean;
    }
    return offset;
}

double HerEmpiricClone::setNewBirth(double mean, double var){
    double offset = 0;
    if (is_mult){
        offset = drawEmpirical(1, var);
        if (offset < 0){
            offset = 0;
        }
        birth_rate = offset * mean;
    }
    else{
        offset = drawEmpirical(0, var);
        if (birth_rate + offset < 0){
            offset = -birth_rate;
        }
        birth_rate = offset + mean;
    }
    return offset;
}

void HerEmpiricClone::reproduce(){
    uniform_real_distribution<double> runif;
    if (runif(*eng) < mut_prob){
        MutationHandler& mut_handle = cell_type->getMutHandler();
        mut_handle.generateMutant(*cell_type, birth_rate, mut_prob);
        removeOneCell();
        double offset = setNewBirth(birth_rate, var);
        HerEmpiricClone *new_node = new HerEmpiricClone(mut_handle.getNewType(), mut_handle.getNewBirthRate(), var, mut_handle.getNewMutProb(), offset, is_mult);
        mut_handle.getNewType().insertClone(*new_node);
        addCells(1);
    }
    else{
        removeOneCell();
        HerEmpiricClone *new_node = new HerEmpiricClone(*cell_type, birth_rate, var, mut_prob, is_mult);
        birth_rate = new_node->getBirthRate();
        addCells(1);
        cell_type->insertClone(*new_node);
    }
}

bool SimpleClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone SimpleClone [type_id] [num_cells] [birth_rate] [mut_rate]
    try{
        cell_count =stoi(parsed_line[0]);
        birth_rate =stod(parsed_line[1]);
        mut_prob =stod(parsed_line[2]);
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>3){
        double death = stod(parsed_line[3]);
        cell_type->setDeathRate(death);
    }
    return checkRep();
}

bool TypeSpecificClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone TypeSpecificClone [type_id] [num_cells] [mean] [var] [mut_rate]
    cell_count = 1;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>4){
        double death = stod(parsed_line[4]);
        cell_type->setDeathRate(death);
    }
    birth_rate = drawFromDist(mean, var);
    return checkRep();
}

bool TypeEmpiricClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone TypeEmpiricClone [type_id] [num_cells] [mean] [var] [mut_rate] [dist_filename]
    cell_count = 1;
    string filename;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
        filename = parsed_line[4];
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>5){
        double death = stod(parsed_line[5]);
        cell_type->setDeathRate(death);
    }
    if (!readDist(filename)){
        return false;
    }
    birth_rate = drawEmpirical(mean, var);
    return checkRep();
}

bool HeritableClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HeritableClone [type_id] [num_cells] [mean] [var] [mut_rate]
    cell_count = 1;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>4){
        dist_type = parsed_line[4];
    }
    if (parsed_line.size()>5){
        double death = stod(parsed_line[5]);
        cell_type->setDeathRate(death);
    }
    setNewBirth(mean, var);
    return checkRep();
}

bool HerResetClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HerResetClone [type_id] [num_cells] [mean] [var] [mut_rate] [num_gen]
    cell_count = 1;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
        num_gen_persist =stoi(parsed_line[4]);
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>5){
        dist_type = parsed_line[5];
    }
    if (parsed_line.size()>6){
        double death = stod(parsed_line[6]);
        cell_type->setDeathRate(death);
    }
    double offset = setNewBirth(mean, var);
    for (int i=0; i<num_gen_persist-1; i++){
        active_diff.push(1);
    }
    active_diff.push(offset);
    return Clone::checkRep() && HerResetClone::checkRep();
}

bool HerResetExpClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HerResetExpClone [type_id] [num_cells] [mean] [var] [mut_rate] [alteration_avg_lifetime] [accum_rate]
    cell_count = 1;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
        time_constant = 1/stod(parsed_line[4]);
        accum_rate =stod(parsed_line[5]);
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>6){
        dist_type = parsed_line[6];
    }
    if (parsed_line.size()>7){
        double death = stod(parsed_line[7]);
        cell_type->setDeathRate(death);
    }
    birth_rate = mean;
    add_alterations();
    return Clone::checkRep() && HerResetExpClone::checkRep();
}

bool HerResetEmpiricClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HerResetEmpiric [type_id] [num_cells] [mean] [var] [mut_rate] [num_gen] [filename]
    cell_count = 1;
    string filename;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
        num_gen_persist =stoi(parsed_line[4]);
        filename = parsed_line[5];
    }
    catch (...){
        return false;
    }
    if (!readDist(filename)){
        return false;
    }
    if (parsed_line.size()>6){
        double death = stod(parsed_line[6]);
        cell_type->setDeathRate(death);
    }
    double offset = setNewBirth(mean, var);
    for (int i=0; i<num_gen_persist-1; i++){
        active_diff.push(1);
    }
    active_diff.push(offset);
    return Clone::checkRep() && HerResetEmpiricClone::checkRep();
}

bool HerEmpiricClone::readLine(vector<string>& parsed_line){
    //full line syntax: Clone HerEmpiricClone [type_id] [num_cells] [mean] [var] [mut_rate] [dist_filename]
    cell_count = 1;
    string filename;
    try{
        mean =stod(parsed_line[1]);
        var =stod(parsed_line[2]);
        mut_prob =stod(parsed_line[3]);
        filename = parsed_line[4];
    }
    catch (...){
        return false;
    }
    if (parsed_line.size()>5){
        double death = stod(parsed_line[5]);
        cell_type->setDeathRate(death);
    }
    if (!readDist(filename)){
        return false;
    }
    birth_rate = drawEmpirical(mean, var);
    return checkRep();
}

bool EmpiricalClone::readDist(string filename){
    if (cell_type->hasDist()){
        return true;
    }
    ifstream infile;
    infile.open(filename);
    if (!infile.is_open()){
        return false;
    }
    string line;
    while (getline(infile, line)){
        std::stringstream ss;
        string tok;
        ss.str(line);
        try{
            while (getline(ss, tok, '\t')){
                cell_type->addDistPoint(stod(tok));
            }
        }
        catch(...){
            return false;
        }
    }
    return cell_type->hasDist();
}

double EmpiricalClone::drawEmpirical(double mean, double var){
    uniform_real_distribution<double> runif;
    int index = floor(runif(*eng) * cell_type->getDistSize());
    return (cell_type->getDistByIndex(index) * sqrt(var)) + mean;
}
