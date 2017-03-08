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

using namespace std;

ThreeTypesMutation::ThreeTypesMutation(double m2, double f1, double f2){
    mu2 = m2;
    fit1 = f1;
    fit2 = f2;
    birth_rate = 0;
    mut_prob = 0;
    new_type = NULL;
}

CellType* MutationHandler::getNewTypeByIndex(int index, CellType& curr_type){
    CList *clone_list = &curr_type.getPopulation();
    if (clone_list->getTypeByIndex(index)){
        return clone_list->getTypeByIndex(index);
    }
    else{
        CellType *new_type = new CellType(index, &curr_type);
        clone_list->insertCellType(*new_type);
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
