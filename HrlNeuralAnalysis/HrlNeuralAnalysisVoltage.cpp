/*
    HRLAnalysis(TM) Software License - Version 1.0 - August 27th, 2013

    Permission is hereby granted, free of charge, to any person or 
    organization obtaining a copy of the software and accompanying 
    documentation covered by this license (the "Software") to use, 
    reproduce, display, distribute, execute, and transmit the 
    Software, and to prepare derivative works of the Software, and 
    to permit third-parties to whom the Software is furnished to do 
    so, all subject to the following:

    The copyright notices in the Software and this entire statement, 
    including the above license grant, this restriction and the 
    following disclaimer, must be included in all copies of the 
    Software, in whole or in part, and all derivative works of the 
    Software, unless such copies or derivative works are solely in 
    the form of machine-executable object code generated by a source 
    language processor.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND 
    NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR 
    ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR 
    OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING 
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
    OTHER DEALINGS IN THE SOFTWARE, INCLUDING BUT NOT LIMITED TO THE 
    COMPATIBILITY OF THIS LICENSE WITH OTHER SOFTWARE LICENSES.
*/
#include <HrlNeuralAnalysisVoltage.h>

using namespace hrlAnalysis;
using namespace std;

// Public Functions

HrlNeuralAnalysisVoltage::HrlNeuralAnalysisVoltage():
		HrlNeuralAnalysis(), voltage_(new VoltageInfo()), numNeurons_(0),
        getSpikes_(0), spikeThreshold_(0)  {}

HrlNeuralAnalysisVoltage::HrlNeuralAnalysisVoltage( 
                    int startTimeIn, int endTimeIn, int startIdxIn,
                    int endIdxIn, std::vector<std::string> fileNames,
                    int numNeuronsIn, bool getSpikesIn, float spikeThresholdIn):
                    
                        HrlNeuralAnalysis(startTimeIn,endTimeIn,startIdxIn,endIdxIn,fileNames),
                        voltage_(new VoltageInfo()), numNeurons_(numNeuronsIn),
                        getSpikes_(getSpikesIn), spikeThreshold_(spikeThresholdIn)  {}

HrlNeuralAnalysisVoltage::~HrlNeuralAnalysisVoltage() {
    clearDataStructures();
}

const VoltageInfoPtr HrlNeuralAnalysisVoltage::voltages() {
    if(!paramsIn_->isDataCompiled)
        buildDataStructures();
    return voltage_;
}

bool HrlNeuralAnalysisVoltage::buildDataStructures() {

    assert(sizeof(int) == 4);
    assert(sizeof(uint) == 4);

    ifstream fpIn;
    bool processFileRetVal(true);
    time_ = 0;
    numCells_ = paramsIn_->endIdx - paramsIn_->startIdx + 1; // number of cells analyzed
    numTimes_ = paramsIn_->endTime - paramsIn_->startTime; // number of time slots analyzed

    float *voltageByTime = new float[numCells_]();
    bool retVal = true;

    // Clear out the existing data structures if needed.
    clearDataStructures();
    if(getSpikes_) {
        cellActivity_->resize(paramsIn_->endIdx - paramsIn_->startIdx + 1);
    }

    voltage_->voltage.resize(numCells_);
    for (int t = 0; t < numCells_; t++)
        voltage_->voltage.at(t).reserve(numTimes_);

    // Loop through the input files and read the voltages.
    BOOST_FOREACH( string file, paramsIn_->fileNames ) {
        // Open the binary file
        fpIn.open(file.c_str(), std::ios::in | std::ios::binary);
        if(!fpIn.fail()) {
            processFileRetVal = processFile(fpIn, voltageByTime);
        } else {
            cerr << "Failed to correctly open file: \n\t" << file << endl;
            retVal = false;
        }
        fpIn.close();

        if(!processFileRetVal) {
            break;
        }
    }
    delete voltageByTime;

    paramsIn_->isDataCompiled = true;
    return retVal;

}

// Private Functions
void HrlNeuralAnalysisVoltage::clearDataStructures() {
    HrlNeuralAnalysis::clearDataStructures();
    BOOST_FOREACH( vector<float> cellVolts, voltage_->voltage) {
        cellVolts.clear();
    }
}

bool HrlNeuralAnalysisVoltage::processFile(ifstream &fpIn, float * tempVoltage) {
    bool retVal = true;

    // Iterate through the file until we're at the start time.
        fpIn.seekg(0, fpIn.end);
        int fEnd = fpIn.tellg();
        fpIn.seekg(0, fpIn.beg);
    while( (paramsIn_->startTime > time_) && fpIn.tellg() < fEnd ) {
            fpIn.seekg(sizeof(float)*(numNeurons_), ios::cur);
        ++time_;
    }

    if(fpIn.good()) {
        // Read in the voltage information.
        while( fpIn.seekg(sizeof(float)*(paramsIn_->startIdx), ios::cur) &&
                fpIn.read((char *) tempVoltage, sizeof(float)*(numCells_)) &&
                fpIn.seekg(sizeof(float)*((numNeurons_ - paramsIn_->endIdx - 1)), ios::cur) ) {

            if(!fpIn.good()) {
                cerr << "Error reading from file." << endl;
                return false;
            }

            for(int i = 0; i < numCells_; i++) {
                voltage_->voltage.at(i).push_back(tempVoltage[i]);
                                
                if(getSpikes_) {
                    if(tempVoltage[i] >= spikeThreshold_) {
                        spikeActivity_->push_back( make_pair(time_, i + paramsIn_->startIdx));
                        cellActivity_->at(i).push_back(time_);
                    }
                }
            }

            ++time_;
            if(!(paramsIn_->startTime <= time_ &&  time_ <= paramsIn_->endTime) ) {
                retVal = false;
                break;
            }
        }
    }
    return retVal;
}

void HrlNeuralAnalysisVoltage::save(std::string filename) {
    std::ofstream ofs(filename.c_str(),std::ofstream::out | std::ofstream::binary);
    boost::archive::binary_oarchive ar(ofs);
    ar & *this;
}

void HrlNeuralAnalysisVoltage::load(std::string filename) {
	std::ifstream ifs(filename.c_str(),std::ifstream::in | std::ifstream::binary);
	boost::archive::binary_iarchive ar(ifs);
	ar & *this;
}