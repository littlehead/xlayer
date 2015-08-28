/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2015 Yuping Dong
 *
 * Author: Yuping Dong <yuping.dong@tufts.edu>
 */

#ifndef XLAYER_MODULE_H
#define XLAYER_MODULE_H

#include "types.hpp"


    
    /*typedef struct InitializedOpts {
    std::vector<float> *stateval;
    std::vector<float> *opteps;
    std::vector<float> *opttau;
    std::vector<float> *optome;
    std::vector<float> *opta1;
    std::vector<uint32_t> *opta2;
    std::vector<uint32_t> *opta3;
    std::vector<uint32_t> *optb1;
} InitializedOpts;*/
    
class CrossLayer {
public:
  static float m_a1[10];
  static float m_s1[10];
  static int m_a2[2];
  static int m_a3[3];
  static int m_b1[4];
  static int m_b2;
  
  static float m_s2[5];
  ValueIterationParameters m_Param;
  float m_ThruWght;
  float m_DelWght;
  float m_CstWght;
  static APPST m_s3[25];
    /*std::vector<float> *m_stateval;
    std::vector<float> *m_opteps;
    std::vector<float> *m_opttau;
    std::vector<float> *m_optome;
    std::vector<float> *m_opta1;
    std::vector<uint32_t> *m_opta2;
    std::vector<uint32_t> *m_opta3;
    std::vector<uint32_t> *m_optb1;
    InitializedOpts m_initopts;*/
  float m_stateval[10][5][25];
  float m_opteps[10][5][25];
  float m_opttau[10][5][25];
  float m_optome[10][5][25];
  float m_opta1[10][5][25];
  uint32_t m_opta2[10][5][25];
  uint32_t m_opta3[10][5][25];
  uint32_t m_optb1[10][5][25];

  CrossLayer ();
  CrossLayer (float PacketTime, int PacketSize, float ThroughputWeight, float DelayWeight, float CostWeight);
  void Optimizer ();
  RETOPT *GetOpt (float PhyState, float MacState, APPST AppState);
};

//}	// namespace ns3

#endif
