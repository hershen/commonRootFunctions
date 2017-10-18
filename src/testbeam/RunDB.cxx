#include "testbeam/RunDB.h"
#include <iostream>
using namespace myFuncs::testbeam;

//-----------------------------------------------------------
// Singleton class to access DB
//-----------------------------------------------------------
RunDB::RunDB()
    : m_DB({{566, RunParams(566, Crystal::CsI_Tl_Belle, -1.0, 58, 0.0, 0.0, 0.085)},
            {567, RunParams(567, Crystal::CsI_Tl_Belle, -1.0, 58, 300, 329.146, 0.085)},
            {568, RunParams(568, Crystal::CsI_Tl_Belle, -1.0, 58, 300, 329.146, 0.085)},
            {570, RunParams(570, Crystal::CsI_Tl_Belle, -1.0, 58, -140, -165.8021138, 0.087)},
            {571, RunParams(571, Crystal::CsI_Tl_Belle, -1.0, 58, -140, -165.8021138, 0.087)},
            {572, RunParams(572, Crystal::CsI_Tl_Belle, -1.0, 58, 300, 329.146, 0.087)},
            {573, RunParams(573, Crystal::CsI_Tl_Belle, -1.0, 58, -100, -124.1836668, 0.087)},
            {574, RunParams(574, Crystal::CsI_Tl_Belle, -1.0, 58, -120, -142.8936974, 0.087)},
            {575, RunParams(575, Crystal::CsI_Tl_Belle, 31.0, 58, -140, -165.8021138, 0.087)},
            {576, RunParams(576, Crystal::CsI_Tl_Belle, 31.0, 58, 0.0, 0.0, 0.087)},
            {577, RunParams(577, Crystal::CsI_Tl_Belle, 31.0, 58, 300, 329.146, 0.087)},
            {578, RunParams(578, Crystal::CsI_Tl_Belle, 6.0, 58, 0.0, 0.0, 0.087)},
            {579, RunParams(579, Crystal::CsI_Tl_Belle, 6.0, 58, -140, -165.8021138, 0.087)},
            {580, RunParams(580, Crystal::CsI_Tl_Belle, 6.0, 58, 300, 329.146, 0.087)},
            {584, RunParams(584, Crystal::CsI_Tl_Belle, 3.0, 58, 0.0, 0.0, 0.087)},
            {585, RunParams(585, Crystal::CsI_Tl_Belle, 3.0, 58, 300, 329.146, 0.087)},
            {586, RunParams(586, Crystal::CsI_Tl_Belle, 3.0, 58, -140, -165.8021138, 0.087)},
            {588, RunParams(588, Crystal::CsI_Tl_Belle, 0.5, 58, 0.0, 0.0, 0.087)},
            {589, RunParams(589, Crystal::CsI_Tl_Belle, 0.5, 58, 300, 329.146, 0.087)},
            {591, RunParams(591, Crystal::CsI_Tl_Belle, 0.5, 58, -140, -165.8021138, 0.087)},
            {592, RunParams(592, Crystal::CsI_Tl_Belle, 3.0, 58, 300, 329.146, 0.087)},
            {593, RunParams(593, Crystal::CsI_Tl_Belle, 6.0, 58, 300, 329.146, 0.087)},
            {594, RunParams(594, Crystal::CsI_Tl_Belle, 31.0, 58, 300, 329.146, 0.087)},
            {595, RunParams(595, Crystal::CsI_Tl_Belle, 31.0, 58, -140, -165.8021138, 0.087)},
            {596, RunParams(596, Crystal::CsI_Tl_Belle, -1.0, 58, 300, 329.146, 0.087)},
            {597, RunParams(597, Crystal::CsI_Tl_Belle, -1.0, 58, -140, -165.8021138, 0.087)},
            {599, RunParams(599, Crystal::CsI_Tl_Belle, -1.0, 58, -120, -142.8936974, 0.087)},
            {600, RunParams(600, Crystal::CsI_Tl_Belle, -1.0, 58, -100, -124.1836668, 0.087)},
            {602, RunParams(602, Crystal::CsI_Tl_Belle, -1.0, 58, 0.0, 0.0, 0.087)},
            {611, RunParams(611, Crystal::CsI_Chinese, -1.0, 458, -120, -142.8936974, 0.082)},
            {612, RunParams(612, Crystal::CsI_Chinese, -1.0, 458, 0.0, 0.0, 0.082)},
            {625, RunParams(625, Crystal::CsI_Chinese, -1.0, 458, -140, -165.8021138, 0.082)},
            {626, RunParams(626, Crystal::CsI_Chinese, -1.0, 458, -140, -165.8021138, 0.082)},
            {627, RunParams(627, Crystal::CsI_Chinese, -1.0, 458, -100, -124.1836668, 0.082)},
            {628, RunParams(628, Crystal::CsI_Chinese, -1.0, 458, 300, 329.146, 0.082)},
            {629, RunParams(629, Crystal::CsI_Chinese, -1.0, 458, 300, 329.146, 0.082)},
            {630, RunParams(630, Crystal::CsI_Chinese, 31.0, 458, 0.0, 0.0, 0.082)},
            {631, RunParams(631, Crystal::CsI_Chinese, 31.0, 458, 300, 329.146, 0.082)},
            {634, RunParams(634, Crystal::CsI_Chinese, 31.0, 458, -140, -165.8021138, 0.082)},
            {635, RunParams(635, Crystal::CsI_Chinese, 6.0, 458, 0.0, 0.0, 0.082)},
            {637, RunParams(637, Crystal::CsI_Chinese, 6.0, 458, -140, -165.8021138, 0.082)},
            {638, RunParams(638, Crystal::CsI_Chinese, 6.0, 458, 300, 329.146, 0.082)},
            {641, RunParams(641, Crystal::CsI_Chinese, 6.0, 458, 300, 329.146, 0.082)},
            {642, RunParams(642, Crystal::CsI_Chinese, 3.0, 458, 0.0, 0.0, 0.082)},
            {643, RunParams(643, Crystal::CsI_Chinese, 3.0, 458, 300, 329.146, 0.082)},
            {646, RunParams(646, Crystal::CsI_Chinese, 3.0, 458, 300, 329.146, 0.082)},
            {650, RunParams(650, Crystal::CsI_Chinese, 3.0, 458, -140, -165.8021138, 0.082)},
            {651, RunParams(651, Crystal::CsI_Chinese, 0.5, 458, 0.0, 0.0, 0.082)},
            {652, RunParams(652, Crystal::CsI_Chinese, 0.5, 458, -140, -165.8021138, 0.082)},
            {653, RunParams(653, Crystal::CsI_Chinese, 0.5, 458, -140, -165.8021138, 0.082)},
            {654, RunParams(654, Crystal::CsI_Chinese, 0.5, 458, 300, 329.146, 0.082)},
            {655, RunParams(655, Crystal::CsI_Chinese, -1.0, 458, 0.0, 0.0, 0.082)},
            {661, RunParams(661, Crystal::CsI_Tl_Babar, -1.0, 458, 0.0, 0.0, 0.0835)},
            {662, RunParams(662, Crystal::CsI_Tl_Babar, -1.0, 458, -100, -124.1836668, 0.0835)},
            {663, RunParams(663, Crystal::CsI_Tl_Babar, -1.0, 458, -120, -142.8936974, 0.0835)},
            {664, RunParams(664, Crystal::CsI_Tl_Babar, -1.0, 458, -120, -142.8936974, 0.0835)},
            {666, RunParams(666, Crystal::CsI_Tl_Babar, -1.0, 458, -140, -165.8021138, 0.0835)},
            {676, RunParams(676, Crystal::CsI_Tl_Babar, -1.0, 458, 300, 329.146, 0.0835)},
            {677, RunParams(677, Crystal::CsI_Tl_Babar, 31.0, 458, -140, -165.8021138, 0.0835)},
            {678, RunParams(678, Crystal::CsI_Tl_Babar, 31.0, 458, 0.0, 0.0, 0.0835)},
            {679, RunParams(679, Crystal::CsI_Tl_Babar, 31.0, 458, 300, 329.146, 0.0835)},
            {680, RunParams(680, Crystal::CsI_Tl_Babar, 6.0, 458, 0.0, 0.0, 0.0835)},
            {681, RunParams(681, Crystal::CsI_Tl_Babar, 6.0, 458, 300, 329.146, 0.0835)},
            {682, RunParams(682, Crystal::CsI_Tl_Babar, 6.0, 458, -140, -165.8021138, 0.0835)},
            {683, RunParams(683, Crystal::CsI_Tl_Babar, 3.0, 458, 0.0, 0.0, 0.0835)},
            {684, RunParams(684, Crystal::CsI_Tl_Babar, 3.0, 458, -140, -165.8021138, 0.0835)},
            {685, RunParams(685, Crystal::CsI_Tl_Babar, 3.0, 458, 300, 329.146, 0.0835)},
            {686, RunParams(686, Crystal::CsI_Tl_Babar, 3.0, 458, 300, 329.146, 0.0835)},
            {688, RunParams(688, Crystal::CsI_Tl_Babar, 0.5, 458, 300, 329.146, 0.0835)},
            {687, RunParams(687, Crystal::CsI_Tl_Babar, 0.5, 458, -140, -165.8021138, 0.0835)},
            {689, RunParams(689, Crystal::CsI_Tl_Babar, -1.0, 458, 0.0, 0.0, 0.0835)},
            {696, RunParams(696, Crystal::CsI_Tl_Babar, -1.0, 458, 0.0, 0.0, 0.0835)},
            {707, RunParams(707, Crystal::CsI_Tl_Babar, -1.0, 458, 300, 329.146, 0.0835)},
            {712, RunParams(712, Crystal::CsI_Ukrainian, -1.0, 458, 300, 329.146, 0.085)},
            {713, RunParams(713, Crystal::CsI_Ukrainian, -1.0, 458, -100, -124.1836668, 0.085)},
            {714, RunParams(714, Crystal::CsI_Ukrainian, -1.0, 458, 0.0, 0.0, 0.085)},
            {715, RunParams(715, Crystal::CsI_Ukrainian, -1.0, 458, -120, -142.8936974, 0.085)},
            {716, RunParams(716, Crystal::CsI_Ukrainian, -1.0, 458, -140, -165.8021138, 0.085)},
            {717, RunParams(717, Crystal::CsI_Ukrainian, 31.0, 458, 300, 329.146, 0.085)},
            {719, RunParams(719, Crystal::CsI_Ukrainian, 31.0, 458, 0.0, 0.0, 0.085)},
            {720, RunParams(720, Crystal::CsI_Ukrainian, 31.0, 458, -140, -165.8021138, 0.085)},
            {722, RunParams(722, Crystal::CsI_Ukrainian, 6.0, 458, 0.0, 0.0, 0.085)},
            {721, RunParams(721, Crystal::CsI_Ukrainian, 6.0, 458, -140, -165.8021138, 0.085)},
            {723, RunParams(723, Crystal::CsI_Ukrainian, 6.0, 458, 300, 329.146, 0.085)},
            {724, RunParams(724, Crystal::CsI_Ukrainian, 3.0, 458, 300, 329.146, 0.085)},
            {725, RunParams(725, Crystal::CsI_Ukrainian, 3.0, 458, 0.0, 0.0, 0.085)},
            {726, RunParams(726, Crystal::CsI_Ukrainian, 3.0, 458, -140, -165.8021138, 0.085)},
            {729, RunParams(729, Crystal::CsI_Ukrainian, 0.5, 458, 0.0, 0.0, 0.085)},
            {730, RunParams(730, Crystal::CsI_Ukrainian, 0.5, 458, 300, 329.146, 0.085)},
            {731, RunParams(731, Crystal::CsI_Ukrainian, 0.5, 458, -140, -165.8021138, 0.085)},
            {732, RunParams(732, Crystal::CsI_Ukrainian, 31.0, 458, 300, 329.146, 0.085)},
            {734, RunParams(734, Crystal::CsI_Ukrainian, -1.0, 458, -140, -165.8021138, 0.085)},
            {735, RunParams(735, Crystal::CsI_Ukrainian, -1.0, 458, -140, -165.8021138, 0.085)},
            {738, RunParams(738, Crystal::CsI_Ukrainian, -1.0, 1240, -140, -165.8021138, 0.085)},
            {737, RunParams(737, Crystal::CsI_Ukrainian, -1.0, 458, -140, -165.8021138, 0.085)},
            {739, RunParams(739, Crystal::CsI_Ukrainian, 31.0, 458, 300, 329.146, 0.085)},
            {743, RunParams(743, Crystal::CsI_Ukrainian, -1.0, 458, 0.0, 0.0, 0.085)},
            {746, RunParams(746, Crystal::CsI_Tl_Belle, -1.0, 58, 300, 329.146, 0.086)},
            {747, RunParams(747, Crystal::CsI_Tl_Belle, 31.0, 58, 300, 329.146, 0.086)},
            {748, RunParams(748, Crystal::CsI_Tl_Belle, 6.0, 58, 300, 329.146, 0.086)},
            {751, RunParams(751, Crystal::CsI_Tl_Belle, 3.0, 58, 300, 329.146, 0.086)},
            {752, RunParams(752, Crystal::CsI_Tl_Belle, 0.5, 58, 300, 329.146, 0.086)},
            {754, RunParams(754, Crystal::CsI_Chinese, -1.0, 458, 0.0, 0.0, 0.092)},
            {755, RunParams(755, Crystal::CsI_Chinese, 3.0, 458, 300, 329.146, 0.092)}} // unordered_map
           )                                                                            // initialization list
{}                                                                                      // constructor

const RunParams &RunDB::operator[](const int runNum) const {
  try {
    return m_DB.at(runNum);
  } catch (std::out_of_range) {
    std::cout << "RunParams::[]: no run with number " << runNum << " in DB" << std::endl;
    throw;
  }
}
