//
//  MutationHandler.cpp
//  evo_sim
//
//  Created by Debra Van Egeren on 2/23/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#include "MutationHandler.h"
#include "main.h"
#include "CList.h"
#include <sstream>
#include <vector>
#include <string>
#include <random>

using namespace std;

MutationHandler::MutationHandler(){
    birth_rate = 0;
    mut_prob = 0;
    has_mutated = false;
    new_type = NULL;
}

ThreeTypesMutation::ThreeTypesMutation(double m2, double f1, double f2) : MutationHandler(){
    mu2 = m2;
    fit1 = f1;
    fit2 = f2;
}

ThreeTypesFlexMutation::ThreeTypesFlexMutation(double m2, double f1, double f2, double pr1) : MutationHandler(){
    mu2 = m2;
    fit1 = f1;
    fit2 = f2;
    p1 = pr1;
}

ManyTypesFlexMutation::ManyTypesFlexMutation(double m2, double f1, double f2, double pr1, int ntypes) : MutationHandler(){
    mu2 = m2;
    fit1 = f1;
    fit2 = f2;
    p1 = pr1;
    num_types  = ntypes;
}

CellType* MutationHandler::getNewTypeByIndex(int index, CellType& curr_type){
    CList *clone_list = &curr_type.getPopulation();
    if (clone_list->getTypeByIndex(index)){
        clone_list->getTypeByIndex(index)->setParent(&curr_type);
        return clone_list->getTypeByIndex(index);
    }
    else{
        CellType *new_type = new CellType(index, &curr_type);
        clone_list->insertCellType(*new_type);
        curr_type.addChild(*new_type);
        return new_type;
    }
}

void ThreeTypesMutation::generateMutant(CellType& type, double b, double mut){
    if (!(type.getIndex() <= 1)){
        throw "bad three types mutating cell type";
    }
    if (type.getIndex() == 1){
        birth_rate = b + fit2;
        mut_prob = 0;
        new_type = getNewTypeByIndex(2, type);
    }
    else if (type.getIndex() == 0){
        birth_rate = b + fit1;
        mut_prob = mu2;
        new_type = getNewTypeByIndex(1, type);
    }
    has_mutated = true;
}

void ThreeTypesFlexMutation::generateMutant(CellType& type, double b, double mut){
    if (!(type.getIndex() <= 1)){
        throw "bad three types mutating cell type";
    }
    if (type.getIndex() == 1){
        birth_rate = fit2;
        mut_prob = 0;
        new_type = getNewTypeByIndex(2, type);
    }
    else if (type.getIndex() == 0){
        uniform_real_distribution<double> runif;
        double which_trans = runif(*eng);
        if (which_trans < p1){
            birth_rate = fit2;
            mut_prob = 0;
            new_type = getNewTypeByIndex(2, type);
        }
        else{
            birth_rate = fit1;
            mut_prob = mu2;
            new_type = getNewTypeByIndex(1, type);
        }
    }
    has_mutated = true;
}

void ManyTypesFlexMutation::generateMutant(CellType& type, double b, double mut){
    if (floor(type.getIndex()/num_types) == 1){
        birth_rate = fit2;
        mut_prob = 0;
        new_type = getNewTypeByIndex(2, type);
    }
    else if (floor(type.getIndex()/num_types) == 0){
        uniform_real_distribution<double> runif;
        double which_trans = runif(*eng);
        if (which_trans < p1){
            birth_rate = fit2;
            mut_prob = 0;
            new_type = getNewTypeByIndex(2*num_types + type.getIndex(), type);
        }
        else{
            birth_rate = fit1;
            mut_prob = mu2;
            new_type = getNewTypeByIndex(num_types + type.getIndex(), type);
        }
    }
    else{
        throw "bad mutating cell type";
    }
    has_mutated = true;
}

void ThreeTypesMultMutation::generateMutant(CellType& type, double b, double mut){
    if (!(type.getIndex() <= 1)){
        throw "bad three types mutating cell type";
    }
    if (type.getIndex() == 1){
        birth_rate = b*fit2/fit1;
        mut_prob = 0;
        new_type = getNewTypeByIndex(2, type);
    }
    else if (type.getIndex() == 0){
        birth_rate = b*fit1;
        mut_prob = mu2;
        new_type = getNewTypeByIndex(1, type);
    }
    has_mutated = true;
}

void NeutralMutation::generateMutant(CellType& type, double b, double mut){
    if (type.getPopulation().noTypesLeft()){
        throw "tried to get new type when no types left";
    }
    else{
        new_type = getNewTypeByIndex(type.getPopulation().getNextType(), type);
        birth_rate = b;
        mut_prob = mut;
    }
    has_mutated = true;
}

bool ThreeTypesMutation::read(std::vector<string>& params){
    
    string pre;
    string post;
    bool isMu2 = false;
    bool isFit1 = false;
    bool isFit2 = false;
    for (int i=0; i<int(params.size()); i++){
        string tok = params[i];
        stringstream ss;
        ss.str(tok);
        getline(ss, pre, ',');
        if (!getline(ss, post)){
            return false;
        }
        if (pre=="mu2"){
            isMu2 = true;
            mu2 =stod(post);
        }
        else if (pre=="fit1"){
            isFit1 = true;
            fit1 =stod(post);
        }
        else if (pre=="fit2"){
            isFit2 = true;
            fit2 =stod(post);
        }
        else{
            return false;
        }
    }
    if (!isMu2 || !isFit1 || !isFit2){
        return false;
    }
    return true;
}

bool ThreeTypesFlexMutation::read(std::vector<string>& params){
    
    string pre;
    string post;
    bool isMu2 = false;
    bool isFit1 = false;
    bool isFit2 = false;
    bool isP1 = false;
    for (int i=0; i<int(params.size()); i++){
        string tok = params[i];
        stringstream ss;
        ss.str(tok);
        getline(ss, pre, ',');
        if (!getline(ss, post)){
            return false;
        }
        if (pre=="mu2"){
            isMu2 = true;
            mu2 =stod(post);
        }
        else if (pre=="fit1"){
            isFit1 = true;
            fit1 =stod(post);
        }
        else if (pre=="fit2"){
            isFit2 = true;
            fit2 =stod(post);
        }
        else if (pre=="p1"){
            isP1 = true;
            p1 =stod(post);
        }
        else{
            return false;
        }
    }
    if (!isMu2 || !isFit1 || !isFit2 || !isP1){
        return false;
    }
    return true;
}

bool ManyTypesFlexMutation::read(std::vector<string>& params){
    
    string pre;
    string post;
    bool isMu2 = false;
    bool isNum = false;
    bool isFit1 = false;
    bool isFit2 = false;
    bool isP1 = false;
    for (int i=0; i<int(params.size()); i++){
        string tok = params[i];
        stringstream ss;
        ss.str(tok);
        getline(ss, pre, ',');
        if (!getline(ss, post)){
            return false;
        }
        if (pre=="mu2"){
            isMu2 = true;
            mu2 =stod(post);
        }
        else if (pre=="fit1"){
            isFit1 = true;
            fit1 =stod(post);
        }
        else if (pre=="num"){
            isNum = true;
            num_types =stod(post);
        }
        else if (pre=="fit2"){
            isFit2 = true;
            fit2 =stod(post);
        }
        else if (pre=="p1"){
            isP1 = true;
            p1 =stod(post);
        }
        else{
            return false;
        }
    }
    if (!isMu2 || !isFit1 || !isFit2 || !isP1 || !isNum){
        return false;
    }
    return true;
}

FixedSitesMutation::FixedSitesMutation() : MutationHandler(){
    fitnesses = NULL;
    adj_mat = NULL;
    max_types = 0;
    is_mult = false;
}

FixedSitesMutation::~FixedSitesMutation(){
    if (fitnesses){
        delete fitnesses;
    }
    if (adj_mat){
        delete adj_mat;
    }
}

bool FixedSitesMutation::read(std::vector<string>& params){
    // syntax: sim_params mut_handler_params [max types] [is mult (integer to bool)] [fitnesses file name] [adj matrix file name]
    
    if (!(params.size() == 4)){
        return false;
    }
    max_types = stoi(params[0]);
    fitnesses = new double[max_types];
    adj_mat = new int[max_types * max_types];
    
    is_mult = stoi(params[1]);
    
    ifstream infile;
    infile.open(params[2]);
    if (!infile.is_open()){
        return false;
    }
    string line;
    int count_type = 0;
    while (getline(infile, line)){
        try{
            fitnesses[count_type] = stod(line);
        }
        catch(...){
            return false;
        }
        count_type ++;
    }
    if (count_type != max_types){
        return false;
    }
    
    ifstream infile2;
    infile2.open(params[3]);
    if (!infile2.is_open()){
        return false;
    }
    count_type = 0;
    while (getline(infile2, line)){
        std::stringstream ss;
        string tok;
        ss.str(line);
        int count_line = 0;
        try{
            while (getline(ss, tok, ',')){
                int to_add = stoi(tok);
                if (to_add > max_types){
                    return false;
                }
                adj_mat[getAdjIndex(count_type, count_line)] = to_add;
                count_line ++;
            }
            for (int i=count_line; i<max_types; i++){
                adj_mat[getAdjIndex(count_type, i)] = -1;
            }
        }
        catch(...){
            return false;
        }
        count_type ++;
    }
    if (count_type > max_types){
        return false;
    }
    for (int i=count_type; i<max_types; i++){
        for (int j=0; j<max_types; j++){
            adj_mat[getAdjIndex(i, j)] = -1;
        }
    }

    return true;
}

int FixedSitesMutation::getAdjIndex(int cell_type, int adj_type){
    return cell_type * max_types + adj_type;
}

void FixedSitesMutation::generateMutant(CellType &type, double b, double mut){
    int count_new = 0;
    int orig_type_id = type.getIndex();
    for (int i=0; i<max_types; i++){
        if (adj_mat[getAdjIndex(orig_type_id, i)] > 0){
            count_new ++;
        }
    }
    if (count_new == 0){
        birth_rate = b;
        mut_prob = 0;
        new_type = &type;
        return;
    }
    uniform_real_distribution<double> runif;
    double which_trans = count_new * runif(*eng);
    int new_type_id = adj_mat[getAdjIndex(orig_type_id, floor(which_trans))];
    new_type = getNewTypeByIndex(new_type_id, type);
    mut_prob = mut;
    if (is_mult){
        birth_rate = b * fitnesses[new_type_id]/fitnesses[orig_type_id];
    }
    else{
        birth_rate = b + fitnesses[new_type_id] - fitnesses[orig_type_id];
    }
    has_mutated = true;
}
