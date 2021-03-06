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
#include <SynchronyInfo.h>

using namespace hrlAnalysis;

CellSynchronyInfo::CellSynchronyInfo(std::vector<int> *spikeTrainIn, int startTime, int endTime) :
        ptr_t_p1_(spikeTrainIn->begin()), ptr_t_p_(spikeTrainIn->begin()),
        ptr_t_f_(spikeTrainIn->begin()), ptr_t_f1_(spikeTrainIn->begin()+1),
        t_p1_(0), t_p_(0), t_f_(0), t_f1_(0),
        dt_p_(0), dt_f_(0), x_isi_(0), x_p_(0), x_f_(0),
        T_(startTime), spikeTrain_(spikeTrainIn), startTime_(startTime), endTime_(endTime)  {};

CellSynchronyInfo::~CellSynchronyInfo() {}

// This has to be called in order.
void CellSynchronyInfo::incrementSynchronyTiming() {
    T_++;
    // Deal with the begining of the spike train when T is before any actual spikes occured.
    if(T_ < *ptr_t_p_ && ptr_t_p_ == spikeTrain_->begin()) {
        t_p1_ = startTime_;
        t_p_ = startTime_;
        t_f_ = *ptr_t_f_;
        t_f1_ = *ptr_t_f1_;
    } else {

        if(ptr_t_p_ != spikeTrain_->end()) {
            if( (t_p_ == startTime_) && (*ptr_t_p_ <= T_) ) {
                t_p_ = *ptr_t_p_;
            } else {
                if( *(ptr_t_p_ + 1) <= T_ && (ptr_t_p_ + 1) != spikeTrain_->end()) {
                    if(ptr_t_p_ != spikeTrain_->begin()) {
                        ptr_t_p1_++;
                    }
                    t_p1_ = *ptr_t_p1_;
                    ptr_t_p_++;
                    t_p_ = *ptr_t_p_;
                }
            }
        }

        if( ptr_t_f_ != spikeTrain_->end() && *(ptr_t_f_) <= T_) {
            ptr_t_f_++;
            if(ptr_t_f_ == spikeTrain_->end()) {
                t_f_ = endTime_;
            } else {
                t_f_ = *ptr_t_f_;
            }

            if(ptr_t_f1_ != spikeTrain_->end()) {
                ptr_t_f1_++;
                if(ptr_t_f1_ == spikeTrain_->end()) {
                    t_f1_ = endTime_;
                } else {
                    t_f1_ = *ptr_t_f1_;
                }

            }

        }

    }

    x_p_ = (double) (T_ - t_p_);
    x_f_ = (double) (t_f_ - T_);
    x_isi_ = (double) (t_f_ - t_p_);

}

// This is not ideal.
void CellSynchronyInfo::calcDeltas(boost::shared_ptr<CellSynchronyInfo> cell2Info) {
    std::vector<int> dt_p_array;
    dt_p_array.push_back(abs(t_p_ - cell2Info->t_p1_));
    dt_p_array.push_back(abs(t_p_ - cell2Info->t_p_));
    dt_p_array.push_back(abs(t_p_ - cell2Info->t_f_));
    dt_p_array.push_back(abs(t_p_ - cell2Info->t_f1_));
    dt_p_ = *std::min_element( dt_p_array.begin(), dt_p_array.end() );

    std::vector<int> dt_f_array;
    dt_f_array.push_back(abs(t_f_ - cell2Info->t_p1_));
    dt_f_array.push_back(abs(t_f_ - cell2Info->t_p_));
    dt_f_array.push_back(abs(t_f_ - cell2Info->t_f_));
    dt_f_array.push_back(abs(t_f_ - cell2Info->t_f1_));
    dt_f_ = *std::min_element( dt_f_array.begin(), dt_f_array.end() );
}

double CellSynchronyInfo::getSn() {
    return ( (double)(dt_p_ * x_f_ + dt_f_ * x_p_) ) / ((double) x_isi_ );
}

















