//
//  OutputWriter.cpp
//  evo_sim
//
//  Created by Debra Van Egeren on 2/21/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#include "OutputWriter.h"
#include "main.h"
#include "Clone.h"
#include "CList.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

FinalOutputWriter::FinalOutputWriter(string ofile){
    ofile_loc = ofile;
    sim_number = 1;
}

DuringOutputWriter::DuringOutputWriter(string ofile, int period){
    ofile_loc = ofile;
    last_written = 0;
    writing_period = period;
}

DuringOutputWriter::DuringOutputWriter(string ofile){
    ofile_loc = ofile;
    last_written = 0;
    writing_period = 0;
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
}

void TypeStructureWriter::beginAction(CList &clone_list){
    string ofile_middle = "sim_"+to_string(sim_number);
    outfile.open(ofile_loc + ofile_middle + ofile_name);
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
    outfile.flush();
    outfile.close();
    sim_number++;
}

CellCountWriter::CellCountWriter(string ofile, int period, int i, int sim):DuringOutputWriter(ofile, period){
    ofile_name = "type_" + to_string(i) + ".oevo";
    index = i;
    sim_number = sim;
}

CellCountWriter::CellCountWriter(string ofile): DuringOutputWriter(ofile){
    sim_number = 1;
}

bool CellCountWriter::readLine(vector<string>& parsed_line){
    if (parsed_line.size() != 2){
        return false;
    }
    try{
        writing_period = stoi(parsed_line[0]);
        index = stoi(parsed_line[1]);
    }
    catch (const invalid_argument){
        return false;
    }
    ofile_name = "type_" + to_string(index) + ".oevo";
    return true;
}

void CellCountWriter::beginAction(CList& clone_list){
    string ofile_middle = "count_sim_"+to_string(sim_number);
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
    ofile_name = "all_types";
}

AllTypesWriter::AllTypesWriter(string ofile): DuringOutputWriter(ofile){
    sim_number = 1;
    ofile_name = "all_types";
}

void AllTypesWriter::beginAction(CList& clone_list){
    int type_index;
    vector<CellType *> root_types = clone_list.getRootTypes();
    for (vector<CellType *>::iterator it = root_types.begin(); it != root_types.end(); ++it){
        type_index = (*it)->getIndex();
        CellCountWriter *new_writer = new CellCountWriter(ofile_loc, writing_period, type_index, sim_number);
        writers.push_back(new_writer);
        new_writer->beginAction(clone_list);
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

bool AllTypesWriter::readLine(vector<string>& parsed_line){
    if (parsed_line.size() != 1){
        return false;
    }
    try{
        writing_period = stoi(parsed_line[0]);
    }
    catch (const invalid_argument){
        return false;
    }
    return true;
}

IsExtinctWriter::IsExtinctWriter(string ofile): FinalOutputWriter(ofile){
    ofile_name = "extinction.oevo";
    outfile.open(ofile_loc+ofile_name);
}

void IsExtinctWriter::finalAction(CList& clone_list){
    outfile << sim_number << ", " << clone_list.isExtinct() << endl;
    sim_number++;
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
    outfile << sim_number << ", " << clone_list.getCurrTime() << endl;
    sim_number++;
}

EndTimeWriter::~EndTimeWriter(){
    outfile.flush();
    outfile.close();
}

IfType2Writer::IfType2Writer(string ofile): FinalOutputWriter(ofile){
    ofile_name = "iftype2.oevo";
    outfile.open(ofile_loc+ofile_name);
}

void IfType2Writer::finalAction(CList& clone_list){
    bool is_2 = clone_list.getTypeByIndex(2);
    is_2 = is_2 && (clone_list.getTypeByIndex(2)->getNumCells() > 0);
    outfile << sim_number << ", " << is_2 << endl;
    sim_number++;
}

IfType2Writer::~IfType2Writer(){
    outfile.flush();
    outfile.close();
}

FitnessDistWriter::FitnessDistWriter(string ofile, int period, int i, int sim):DuringOutputWriter(ofile, period){
    ofile_name = "type_" + to_string(i) + ".oevo";
    index = i;
    sim_number = sim;
}

FitnessDistWriter::FitnessDistWriter(string ofile): DuringOutputWriter(ofile){
    sim_number = 1;
}

bool FitnessDistWriter::readLine(vector<string>& parsed_line){
    if (parsed_line.size() != 2){
        return false;
    }
    try{
        writing_period = stoi(parsed_line[0]);
        index = stoi(parsed_line[1]);
    }
    catch (const invalid_argument){
        return false;
    }
    ofile_name = "type_" + to_string(index) + ".oevo";
    return true;
}

void FitnessDistWriter::beginAction(CList& clone_list){
    string ofile_middle = "fit_sim_"+to_string(sim_number);
    outfile.open(ofile_loc + ofile_middle + ofile_name);
    outfile << "data for cell type " << index << " sim number " << sim_number << endl;
    outfile << clone_list.getCurrTime();
    write_dist(outfile, clone_list);
    outfile << endl;
    
}

FitnessDistWriter::~FitnessDistWriter(){
    outfile.flush();
    outfile.close();
}

void FitnessDistWriter::write_dist(ofstream& outfile, CList& clone_list){
    Clone *curr_clone = (clone_list.getTypeByIndex(index)->getRoot());
    while (curr_clone){
        int num_cells = curr_clone->getCellCount();
        for (int i=0; i<num_cells; i++){
            outfile << ", " << curr_clone->getBirthRate();
        }
        curr_clone = &(curr_clone->getNextWithinType());
    }
}

void FitnessDistWriter::duringSimAction(CList& clone_list){
    if (shouldWrite(clone_list) && clone_list.getTypeByIndex(index)->getNumCells() > 0){
        outfile << clone_list.getCurrTime();
        write_dist(outfile, clone_list);
        outfile << endl;
    }
}

void FitnessDistWriter::finalAction(CList& clone_list){
    outfile << clone_list.getCurrTime();
    write_dist(outfile, clone_list);
    outfile << endl;
    outfile.flush();
    outfile.close();
    sim_number++;
}
