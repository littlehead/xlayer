/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2015 Yuping Dong
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * Authors: Yuping Dong  <yuping.dong@tufts.edu>
 *          
 */


///
/// \brief Implementation of cross-layer design.
///
/// This is the main file of this software.
///


#include <iostream>
#include <fstream>
#include <new>
#include <cmath>
#include <ctime>
#include "xlayer-modulel.hpp"
#include "types.hpp"
#include "basefuncs.hpp"
#include "qosfuncs.hpp"
#include "iterfuncs.hpp"

using namespace std;

/********** Useful macros **********/

///
/// \brief Gets the delay between a given time and the current time.
///
/// If given time is previous to the current one, then this macro returns
/// a number close to 0. This is used for scheduling events at a certain moment.
///
//#define DELAY(time) (((time) < (Simulator::Now ())) ? Seconds (0.000001) : \
                     (time - Simulator::Now () + Seconds (0.000001))


//namespace ns3 {

//NS_LOG_COMPONENT_DEFINE ("CrossLayer");
  

float CrossLayer::m_a1[10];
float CrossLayer::m_s1[10];
int CrossLayer::m_a2[2];
int CrossLayer::m_a3[3];
int CrossLayer::m_b1[4];
int CrossLayer::m_b2;
float CrossLayer::m_s2[5];
APPST CrossLayer::m_s3[25];

CrossLayer::CrossLayer (float PacketTime, int PacketSize, float ThroughputWeight, float DelayWeight, float CostWeight)
{
    int i, j;
    
    
    for (i = 0; i < 10; i++)
    {
        m_a1[i] = 16 + 2 * i;
        m_s1[i] = 8 + 4 * i + 101;
    }
    
    m_a2[0] = 0;
    m_a2[1] = 1;
    
    for (i = 0; i < 3; i++)
    {
        m_a3[i] = m_b1[i] = i + 1;
    }
    
    m_b1[3] = 4;
    m_b2 = 5;
    
    for (i = 0; i < 5; i++)
        //m_s2[i] = 0.2 + 0.2 * i;
        m_s2[i] = 1;
    
    m_Param.lamda1 = 1;
    m_Param.lamda2 = 1;
    m_Param.lamda3 = 1;
    m_Param.lamdag = 0.1;
    m_Param.gama = 0.9;
    m_Param.T = 0.02;
    m_Param.Tp = 0.0008;
    m_Param.PktSize = 4;
    m_Param.DopplerFreq = 50;
 
  m_ThruWght = ThroughputWeight;
  m_DelWght = DelayWeight;
  m_CstWght = CostWeight;
    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < 5; j++)
        {
            m_s3[i * 5 + j].x = i;
            m_s3[i * 5 + j].y = j;
        }
    }

    Optimizer ();
    
}


    
void
CrossLayer::Optimizer ()
{
    int i, j, k, m, x;
    int PktSize;
    float next_val[10][5][25];
    QOS *Z1 = NULL;
    QOS *Z3 = NULL;
    QOS *Z_tmp = NULL;
    RETLY1 *stval_tmp;
    clock_t begin, end;
    double proc_time = 0.0;
    double ave_proc_time;
    float diff, tmp_diff, ave_tmp_diff, act_diff;
    float thruwght, delwght, cstwght, Tp;
    //ofstream qosfile0, qosfile1, qosfile2;
    ofstream stvalfile;
    int test;
    
    thruwght = m_ThruWght;
    delwght = m_DelWght;
    cstwght = m_CstWght;
    PktSize = m_Param.PktSize;
    Tp = m_Param.Tp;
    
    /*std::cout << "possible states S3s are " << std::endl;

    for (i = 0; i < 25; i++)
        std::cout << "{" << CrossLayer::m_s3[i].x << ", " << CrossLayer::m_s3[i].y << "} " ;*/
    
    /*std::cout << "ThroughputWeight is " << thruwght << ", DelayWeight is " << delwght << ", CostWeight is " << cstwght << std::endl;
    std::cout << "PacketSize is " << PktSize << ", PacketTime is " << Tp << std::endl;*/
    
    x = 0;
    
    for (k = 0; k < 10; k++)
    {
        for (j = 0; j < 5; j++)
        {
            for (m = 0; m < 25; m++)
                next_val[k][j][m] = 0;
        }
    }
    
    diff = 10;
    
    stvalfile.open("stval1.txt");
    //std::cout << "epsilon = " << epsilon << std::endl;
    
    while (diff > epsilon)
    {
    /*for (test = 0; test < 3; test++)
    {*/
        /*if (test == 0)
            qosfile0.open ("z10.txt");
        else if (test == 1)
            qosfile1.open ("z11.txt");
        else
            qosfile2.open("z12.txt");*/
        
        tmp_diff = 0;
        
        for (k = 0; k < 10; k++)        //iterate on s1
        {
            //QOS calculation
            for (i = 0; i < 4; i++)   //for each b1, calculate Z1(b1)
                calcqosphy(m_s1[k], m_b1[i], Tp, PktSize, &Z1);
            
            /*Z_tmp = Z1;
            while (Z_tmp != NULL)
            {
                printf("PHY layer QOS pkt loss rate = %f, tx time = %f, tx cost = %f.\n", Z_tmp->x, Z_tmp->y, Z_tmp->z);
                Z_tmp = Z_tmp->next;
            }*/
            
            for (j = 0; j < 5; j++)      //iterate on s2
            {
                /*Z_tmp = Z1;
                 while (Z_tmp != NULL)
                 {
                 printf("PHY layer QOS pkt loss rate = %f, tx time = %f, tx cost = %f.\n", Z_tmp->x, Z_tmp->y, Z_tmp->z);
                 Z_tmp = Z_tmp->next;
                 }*/
                
                //printf("calculating Z3...\n");
                Z3 = optqosfrontier(Z1, m_s2[j], m_b2, thruwght, delwght, cstwght);
                
                /*if (test == 0)
                {
                    Z_tmp = Z3;
                    while (Z_tmp != NULL) {
                        qosfile0 << "s1 = " << m_s1[k] << ",s2 = " << m_s2[j] << ", Z1 = (" << Z_tmp->x << ", " << Z_tmp->y << ", " << Z_tmp->z << ")\n";
                        Z_tmp = Z_tmp->next;
                    }
                }
                else if (test == 1)
                {
                    Z_tmp = Z3;
                    while (Z_tmp != NULL) {
                        qosfile1 << "s1 = " << m_s1[k] << ",s2 = " << m_s2[j] << ", Z1 = (" << Z_tmp->x << ", " << Z_tmp->y << ", " << Z_tmp->z << ")\n";
                        Z_tmp = Z_tmp->next;
                    }
                }
                else
                {
                    Z_tmp = Z3;
                    while (Z_tmp != NULL) {
                        qosfile2 << "s1 = " << m_s1[k] << ",s2 = " << m_s2[j] << ", Z1 = (" << Z_tmp->x << ", " << Z_tmp->y << ", " << Z_tmp->z << ")\n";
                        Z_tmp = Z_tmp->next;
                    }
                }*/
                
                /*Z_tmp = Z3;
                 while (Z_tmp != NULL)
                 {
                 printf("APP layer QOS pkt loss rate = %f, tx time = %f, tx cost = %f, b1 = %d\n", Z_tmp->x, Z_tmp->y, Z_tmp->z, Z_tmp->b1);
                 //printf("average transmit time per packet is %f\n", Z_tmp->y/(1-Z_tmp->x));
                 Z_tmp = Z_tmp->next;
                 }*/
                
                for (m = 0; m < 25; m++)  //iterate on s3
                {
                    //value iteration
                    
                    begin = clock();
                    stval_tmp = stvalfunc(m_s1[k], m_s2[j], Z3, m_s3[m], next_val, thruwght, delwght, cstwght, m_Param);
                    end = clock();
                    proc_time += (double)(end - begin)/CLOCKS_PER_SEC;
                    m_stateval[k][j][m] = stval_tmp->x;
                    m_opteps[k][j][m] = stval_tmp->zx;
                    m_opttau[k][j][m] = stval_tmp->zy;
                    m_optome[k][j][m] = stval_tmp->zz;
                    m_opta1[k][j][m] = stval_tmp->a1;
                    m_opta2[k][j][m] = stval_tmp->a2;
                    m_opta3[k][j][m] = stval_tmp->a3;
                    m_optb1[k][j][m] = stval_tmp->b1;
                    act_diff = m_stateval[k][j][m] - next_val[k][j][m];
                    tmp_diff += fabsf (act_diff);
                    
                }
                deallocmem(Z3);
                Z3 = NULL;
            }
            
            deallocmem(Z1);
            Z1 = NULL;
            
            
        }
        
        for (k = 0; k < 10; k++)        //iterate on s1
        {
            for (j = 0; j < 5; j++)      //iterate on s2
            {
                for (m = 0; m < 25; m++)  //iterate on s3
                    next_val[k][j][m] = m_stateval[k][j][m];
            }
        }
        
        //printf("difference of two iteration is %f\n", tmp_diff);
        
        x++;
        std::cout << x << "th iteration." << std::endl;
        
        ave_tmp_diff = tmp_diff / 1250;
        
        std::cout << "stval difference is " << ave_tmp_diff << "after this iteration." << std::endl;
        
        if(ave_tmp_diff < diff)
            diff = ave_tmp_diff;
        
        
        
    }
    //}
    
    ave_proc_time = proc_time / (x * 10 * 5 * 25);
    //std::cout << "average processing time for each stval calc is " << ave_proc_time << "s" << std::endl;
    for (k = 0; k < 10; k++)        //iterate on s1
    {
        for (j = 0; j < 5; j++)      //iterate on s2
        {
            for (m = 0; m < 25; m++)  //iterate on s3
                stvalfile << m_stateval[k][j][m] << " ";
        }
        stvalfile << "\n";
    }

    
    stvalfile.close();
        /*if (test == 0)
            qosfile0.close();
        else if (test == 1)
            qosfile1.close();
        else
            qosfile2.close();*/
}

RETOPT
*CrossLayer::GetOpt (float PhyState, float MacState, APPST AppState)
{
  RETOPT *optresults = new RETOPT;
  int r, s, t, k, j, m;

  r = s = t = k = j = m = 0;

  while (m_s1){
    if (PhyState == m_s1[k]){
      r = k;
      break;
    } k++;
  }

  while (m_s2){
    if (MacState == m_s2[j]){
      s = j;
      break;
    } j++;
  }

  while (m_s3){
    if ((AppState.x == m_s3[m].x) && (AppState.y == m_s3[m].y)){
      t = m;
      break;
    } m++;
  }

  optresults->value = m_stateval[r][s][t];
  optresults->extaction1 = m_opta1[r][s][t];
  optresults->extaction2 = m_opta2[r][s][t];
  optresults->extaction3 = m_opta3[r][s][t];
  optresults->intaction1 = m_optb1[r][s][t];
  optresults->qosx = m_opteps[r][s][t];
  optresults->qosy = m_opttau[r][s][t];
  optresults->qosz = m_optome[r][s][t];

  return optresults;
}

//} // namespace ns3


