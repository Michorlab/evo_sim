//
//  OutputWriter.cpp
//  evo_sim
//
//  Created by Debra Van Egeren on 2/21/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#include "OutputWriter.h"
#include "main.h"
#include "CList.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

FinalOutputWriter::FinalOutputWriter(string ofile){
    ofile_loc = ofile;
}

DuringOutputWriter::DuringOutputWriter(string ofile, int period){
    ofile_loc = ofile;
    last_written = 0;
    writing_period = period;
}

bool DuringOutputWriter::shouldWrite(CList& clone_list){
    if (writing_period == 0){
        return true;
    }
    int floored_time = floor(clone_list.getCurrTime());
    if ((floored_time % writing_period == 0) && (floored_time != last_written)){
        last_written = floored_time;
        return true;
    }
    return false;
}

TypeStructureWriter::TypeStructureWriter(string ofile):FinalOutputWriter(ofile){
    ofile_name = "type_tree.oevo";
    outfile.open(ofile_loc+ofile_name);
}

TypeStructureWriter::~TypeStructureWriter(){
    outfile.flush();
    outfile.close();
}

void TypeStructureWriter::finalAction(CList& clone_list){
    std::vector<CellType *> roots = clone_list.getRootTypes();
    for (int i=0; i<roots.size(); i++){
        clone_list.walkTypesAndWrite(outfile, *roots[i]);
    }
}

CellCountWriter::CellCountWriter(string ofile, int period, int i, int sim):DuringOutputWriter(ofile, period){
    ofile_name = "type_" + to_string(i) + ".oevo";
    index = i;
    sim_number = sim;
}

void CellCountWriter::beginAction(CList& clone_list){
    string ofile_middle = "sim_"+to_string(sim_number);
    outfile.open(ofile_loc + ofile_middle + ofile_name);
    outfile << "data for cell type " << index << " sim number " << sim_number << endl;
    outfile << clone_list.getCurrTime() << ", " << clone_list.getTypeByIndex(index)->getNumCells() << endl;
}

CellCountWriter::~CellCountWriter(){
    outfile.flush();
    outfile.close();
}

void CellCountWriter::duringSimAction(CList& clone_list){
    if (shouldWrite(clone_list) && clone_list.getTypeByIndex(index)->getNumCells() > 0){
        outfile << clone_list.getCurrTime() << ", " << clone_list.getTypeByIndex(index)->getNumCells() << endl;
    }
}

void CellCountWriter::finalAction(CList& clone_list){
    outfile << clone_list.getCurrTime() << ", " << clone_list.getTypeByIndex(index)->getNumCells() << endl;
    outfile.flush();
    outfile.close();
    sim_number++;
}

AllTypesWriter::AllTypesWriter(string ofile, int period): DuringOutputWriter(ofile, period){
    sim_number = 1;
}

void AllTypesWriter::beginAction(CList& clone_list){
    int type_index;
    vector<CellType *> root_types = clone_list.getRootTypes();
    for (vector<CellType *>::iterator it = root_types.begin(); it != root_types.end(); ++it){
        type_index = (*it)->getIndex();
        CellCountWriter *new_writer = new CellCountWriter(ofile_loc, writing_period, type_index, sim_number);
        writers.push_back(new_writer);
    }
}

void AllTypesWriter::duringSimAction(CList& clone_list){
    int new_index = clone_list.newestType();
    bool made_new_type = true;
    for (vector<CellCountWriter *>::iterator it = writers.begin(); it != writers.end(); ++it){
        (*it)->duringSimAction(clone_list);
        made_new_type = made_new_type && !(new_index == (*it)->getTypeIndex());
    }
    if (made_new_type){
        CellCountWriter *new_writer = new CellCountWriter(ofile_loc, writing_period, new_index, sim_number);
        writers.push_back(new_writer);
        new_writer->beginAction(clone_list);
    }
}

void AllTypesWriter::finalAction(CList& clone_list){
    for (vector<CellCountWriter *>::iterator it = writers.begin(); it != writers.end(); ++it){
        (*it)->finalAction(clone_list);
    }
    writers.clear();
    sim_number++;
}

IsExtinctWriter::IsExtinctWriter(string ofile): FinalOutputWriter(ofile){
    ofile_name = "extinction.oevo";
    outfile.open(ofile_loc+ofile_name);
}

void IsExtinctWriter::finalAction(CList& clone_list){
    outfile << clone_list.isExtinct() << endl;
}

IsExtinctWriter::~IsExtinctWriter(){
    outfile.flush();
    outfile.close();
}

EndTimeWriter::EndTimeWriter(string ofile): FinalOutputWriter(ofile){
    ofile_name = "end_time.oevo";
    outfile.open(ofile_loc+ofile_name);
}

void EndTimeWriter::finalAction(CList& clone_list){
    outfile << clone_list.getCurrTime() << endl;
}

EndTimeWriter::~EndTimeWriter(){
    outfile.flush();
    outfile.close();
}

IfType2Writer::IfType2Writer(string ofile): FinalOutputWriter(ofile){
    ofile_name = "type2.oevo";
    outfile.open(ofile_loc+ofile_name);
}

void IfType2Writer::finalAction(CList& clone_list){
    outfile << (clone_list.getTypeByIndex(2)->getNumCells() > 0) << endl;
}

IfType2Writer::~IfType2Writer(){
    outfile.flush();
    outfile.close();
}
