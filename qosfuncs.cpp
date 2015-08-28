//
//  qosfuncs.c
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>
#include "xlayer-modulel.hpp"
#include "types.hpp"
#include "basefuncs.hpp"
#include "qosfuncs.hpp"

using namespace std;

QOS *calcupperqos(QOS *QOS_lower, float State, int ActionB)
{
    QOS *QOS_up = NULL;
    QOS_up = new QOS;
    QOS_up->x = (float)pow((double)QOS_lower->x, (double)(ActionB + 1));
    QOS_up->y = (1 - (float)pow((double)QOS_lower->x, (double)ActionB)) * QOS_lower->y / ((1 - QOS_lower->x) * State);
    QOS_up->z = (1 - (float)pow((double)QOS_lower->x, (double)ActionB)) * QOS_lower->z / (1 - QOS_lower->x);
    QOS_up->b1 = QOS_lower->b1;
    QOS_up->next = NULL;
    
    return QOS_up;
}

QOS *calcphyqos(float PhyState, int ModLevel, float PktTime, int PktSize)
{
    QOS *QOS_phy = NULL;
    float BER;
    float test;
    int mry;
    
    QOS_phy = new QOS;
    //BER = erfcf(283.5 * PhyState * sin(pi/powf((float)2, (float)ModLevel)));
    test = sqrtf(PhyState) * sin(pi/powf((float)2, (float)ModLevel));
    mry = pow(2, ModLevel);
    BER = (1 / (float)ModLevel) * erfcf(sqrtf(PhyState) * sin(pi/powf((float)2, (float)mry)));
    //printf("BER = %f\n", BER);
    QOS_phy->x = (float)(1 - powf((1 - BER), (float)PktSize));
    QOS_phy->y = (float)(PktTime / ModLevel);
    QOS_phy->z = powf((float)ModLevel/4, 1/3)/500;
    QOS_phy->b1 = ModLevel;
    QOS_phy->next = NULL;
    /*if(ModLevel == 1)
    {
        printf("s1 = %f\n", PhyState);
        printf("erfc(%f) = %f, BER = %f\n", test, erfcf(test), BER);
        printf("PHY layer QOS packet loss rate = %f, tx time = %f, tx cost = %f\n", QOS_phy->x, QOS_phy->y, QOS_phy->z);
    }*/
    return QOS_phy;
}

QOS *optqosfrontier(QOS *QOS_lower, float State, int ActionB, float ThroughputWeight, float DelayWeight, float CostWeight)
{
    int flag = 0;
    QOS *QOS_upper = NULL;
    QOS *QOS_intermediate = NULL;
    QOS *QOS_traverse = NULL;
    QOS *QOS_temp = NULL;
    float thruwght, delwght, cstwght;
    
    thruwght = ThroughputWeight;
    delwght = DelayWeight;
    cstwght = CostWeight;
    //int i = 0;
    //printf("calculating qos frontier...\n");
    while(QOS_lower != NULL)
    {
        //printf("calculating upper qos...\n");
        flag = 0;
        QOS_intermediate = calcupperqos(QOS_lower, State, ActionB);
        
        //printf("QOS_intermediate.x = %f, QOS_intermediate.y = %f, QOS_intermediate.z = %f, QOS_intermediate.b1 = %d\n", QOS_intermediate->x, QOS_intermediate->y, QOS_intermediate->z, QOS_intermediate->b1);
        
        QOS_traverse = QOS_upper;
        QOS_temp = QOS_upper;
        
        while(QOS_traverse != NULL)
        {
            /*if((fabsf(thruwght - delwght) + fabsf(cstwght - delwght)) > 0.01)
            {
                if(((QOS_intermediate->x - QOS_traverse->x) * thruwght / QOS_traverse->x + (QOS_intermediate->y - QOS_traverse->y) * delwght / QOS_traverse->y + (QOS_intermediate->z - QOS_traverse->z) * cstwght / QOS_traverse->z > 0) || ((QOS_intermediate->y / (1 - QOS_intermediate->x)) > 0.005))
                {
                    flag = 1;
                    break;
                }
            }
            else
            {*/
                if((QOS_traverse->x <= QOS_intermediate->x)&&(QOS_traverse->y <= QOS_intermediate->y)&&(QOS_traverse->z <= QOS_intermediate->z))
                {
                    flag = 1;
                    break;
                }
            //}
            QOS_traverse = QOS_traverse->next;
        }
        //printf("adding qos upper...\n");
        
        if(flag == 0)
        {
            add_to_list(QOS_intermediate->x, QOS_intermediate->y, QOS_intermediate->z, QOS_intermediate->b1, &QOS_upper);
            //printf("APP opt QOS.x = %f, QOS.y = %f, QOS.z = %f, b1 = %d.\n", QOS_intermediate->x, QOS_intermediate->y, QOS_intermediate->z, QOS_intermediate->b1);
        }
        /*        i++;
         while(QOS_temp != NULL)
         {
         printf("iteration # %d: QOS_temp.x = %f, QOS_temp.y = %f, QOS_temp.z = %f\n", i, QOS_temp->x, QOS_temp->y, QOS_temp->z);
         QOS_temp = QOS_temp->next;
         }*/
        //printf("calculating next qos upper...\n");
        QOS_lower = QOS_lower->next;
        free(QOS_intermediate);
    }
    
    return QOS_upper;
}



void calcqosphy(float PhyState, int ModLevel, float PktTime, int PktSize, QOS **QOSPhy)
{
    int flag = 0;
    //QOS *QOSPhy = NULL;
    QOS *QOS_traverse = NULL;
    QOS *QOS_intermediate = NULL;
    
    
    flag = 0;
    QOS_intermediate = calcphyqos(PhyState, ModLevel, PktTime, PktSize);
    //printf("PHY layer QOS loss rate = %f, tx time = %f, tx cost = %f\n", QOS_intermediate->x, QOS_intermediate->y, QOS_intermediate->z);
            
    QOS_traverse = *QOSPhy;
            
    while(QOS_traverse != NULL)
    {
        if((QOS_traverse->x <= QOS_intermediate->x)&&(QOS_traverse->y <= QOS_intermediate->y)&&(QOS_traverse->z <= QOS_intermediate->z))
        {
            flag = 1;
            break;
        }
        QOS_traverse = QOS_traverse->next;
    }
            
    if(flag == 0)
        add_to_list(QOS_intermediate->x, QOS_intermediate->y, QOS_intermediate->z, QOS_intermediate->b1, QOSPhy);
            
    free(QOS_intermediate);
    
    //return QOSPhy;
}
