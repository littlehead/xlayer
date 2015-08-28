//
//  iterfuncs.c
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "xlayer-modulel.hpp"
#include "types.hpp"
#include "basefuncs.hpp"
#include "iterfuncs.hpp"
//#include <functional>


using namespace std;

float db2rt(float state)
{
    return (powf(10, (state / 10)));
}

float omega(float state, float ave_SINR)
{
    return (expf(-db2rt(state)/db2rt(ave_SINR)) - expf(-(db2rt(state + 2))/db2rt(ave_SINR)));
}

float xing(float state, float ave_SINR, float DopplerFreq)
{
    return (sqrtf(2 * pi * db2rt(state) / db2rt(ave_SINR)) * DopplerFreq * expf(-db2rt(state)/db2rt(ave_SINR)));
}

RETLY3 *stval2(float next_s1, float next_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param)
{
    float sum_val12 = 0.0;
    float thrupt, delay;
    float Rin;
    float prob3;
    int n3, nLost, nSent, precalcnexts3x, totavailpkts;
    int i, j, k, m;
    APPST next_s3;
    //QOS *Z_opt = NULL;
    RETLY3 *retstruct3 = new RETLY3;
    //RETLY1 *val_12 = (RETLY1 *)malloc(sizeof(RETLY1));
    float thruwght, delwght, cstwght;
    float test;
    
    thruwght = ThroughputWeight;
    delwght = DelayWeight;
    cstwght = CostWeight;
    
    Rin = 0.0;
    retstruct3->x = 0.0;
    
    j = 0;
    k = 0;
    
    //printf("next_s1 = %f, next_s2 = %f.\n",next_s1, next_s2);

    //find next s1, s2
    while (fabs(CrossLayer::m_s1[j] - next_s1) > 0.05)
        j++;
    
    while (fabs(CrossLayer::m_s2[k] - next_s2) > 0.05)
        k++;
    
    /*if ((next_s1 == 4) && (next_s2 == 1))
        cout << "j = " << j << ", k = " << k << ",next_s1 = " << next_s1 << ", next_s2 = " << next_s2 << ", s1[j] = " << CrossLayer::m_s1[j] << ", s2[k] = " << CrossLayer::m_s2[k] << endl;*/
    
    //printf("j = %d, k = %d\n", j, k);
    while(Z3 != NULL)
    {
        n3 = (1 - Z3->x) * Param.T / Z3->y;
        
        /*if ((j == 9) && (k == 4))
            cout << "n3 = " << n3 << endl;*/
        nLost = curr_s3.x - n3;
        //Rin = n3 - lamdag * max(nLost, 0);
        
        //modify Rin here
        thrupt = (n3 - Param.lamdag * max(nLost, 0));  //normalize throughput
        //totavailpkts = curr_s3.x + curr_s3.y;
        //thrupt = min(n3, totavailpkts);
        delay = Z3->y/(1-Z3->x);         //normalize delay for each packet
        //Rin = thruwght * thrupt + delwght * (1/delay);  //add delay to Rin
        Rin = thrupt + 0.0001 * (1/delay);
        //Rin = thrupt;
        
        /*if ((j == 1) && (k == 1))
            cout << "n3 = " << n3 << ", Rin = " << Rin << ".\n";*/
        /*if(Rin < 0)
        {
            printf("s1[%d] = %f, s2[%d] = %f\n", j, s1[j], k, s2[k]);
            printf("plr = %f, tx time = %f and tx cost = %f\n", Z3->x, Z3->y, Z3->z);
            
        }*/
        //cout << "Rin = " << Rin << "with Z3 = { " << Z3->x << ", " << Z3->y << ", " << Z3->z << "\n";
        //printf("Rin = %f with plr = %f, tx time = %f and tx cost = %f\n", Rin, Z3->x, Z3->y, Z3->z);
        for(i = 0; i < 3; i++)      //iterate on a3
        {
            //sum_val12 = Rin - cstwght * lamda3 * a3[i];
            sum_val12 = Rin - Param.lamda3 * CrossLayer::m_a3[i];
            //cout << "sum_val12 is " << sum_val12 << " with Z3 = { " << Z3->x << ", " << Z3->y << ", " << Z3->z << ", " << Z3->b1 << "} and a3 = " << CrossLayer::m_a3[i] << ".\n";
            //printf("sum_val12 is %f with Z3 = {%f, %f, %f, %d} and a3 = %d\n", sum_val12, Z3->x, Z3->y, Z3->z, Z3->b1, a3[i]);
            /*if ((j == 9) && (k == 4))
                cout << "sum_val12 = " << sum_val12 << ", a3 = " << CrossLayer::m_a3[i] << ", Rin = " << Rin << ".\n";*/
           
            
///////*****************************************************************************************************///////////////
///////*****************  Problem is here!! when curr_s3.y = 0, test is larger than 0  *********************///////////////
            
            nSent = n3 - curr_s3.x;
            precalcnexts3x = curr_s3.y - max(nSent, 0);
            //next_s3.x = precalcnexts3x;
            if(precalcnexts3x < 0)
                next_s3.x = 0;
            else
                next_s3.x = precalcnexts3x;
            
            
            //if ((j == 9) && (k == 4))
                //cout << "next_s3.x = " << next_s3.x << "\n";
            //cout << "test = " << test << "\n";
            
            /*if ((j == 9) && (k == 4))
                cout << "test = " << test << "curr_s3.y = " << curr_s3.y << "\n";*/
///////**************************************************************************************************//////////////////
            //printf("next_s3.x = %d\n", next_s3.x);
            //printf("s3[0] = %d %d\n", s3[0].x, s3[0].y);
            
            for(next_s3.y = 0; next_s3.y < 5; next_s3.y++)
            {
                m = 0;
                //printf("next_s3.x = %d, next_s3.y = %d\n", next_s3.x, next_s3.y);
                while ((CrossLayer::m_s3[m].x != next_s3.x)||(CrossLayer::m_s3[m].y != next_s3.y))
                {
                    m++;
                    if(m == 25)
                        break;
                }
        
                //if (i == 0)
                  //  cout << "m = " << m << "\n";
                if (m == 25)
                    prob3 = 0;
                else
                    prob3 = (float)pow((double)CrossLayer::m_a3[i], (double)next_s3.y) * exp((double)(-CrossLayer::m_a3[i])) / fact(next_s3.y);
                
                test = exp((double)(-CrossLayer::m_a3[i]));
                //if (i == 0)
                  //  cout << "exponential is " << test << "\n";
                    //cout << "prob3 = " << prob3 << "\n";
                sum_val12 += Param.gama * prob3 * NEXT_VAL[j][k][m];
                /*if((j < 2)&&(k < 2)&&(m < 2))
                    printf("Next_val[%d][%d][%d] = %f ", j, k, m, NEXT_VAL[j][k][m]);*/
                /*if ((j == 9) && (k == 4))
                    cout << "prob3 is " << prob3 << "\n";*/
                
                //val_12 = stvalfunc(next_s1, next_s2, next_s3, Z3, NEXT_VAL);
                //sum_val12 += gama * prob3 * val_12->x;
            }
            //cout << "sum_val12 is " << sum_val12 << " with Z3 = { " << Z3->x << ", " << Z3->y << ", " << Z3->z << ", " << Z3->b1 << "} and a3 = " << CrossLayer::m_a3[i] << ".\n";
            //printf("sum_val12 is %f with Z3 = {%f, %f, %f, %d} and a3 = %d\n", sum_val12, Z3->x, Z3->y, Z3->z, Z3->b1, a3[i]);
            
            if (sum_val12 > retstruct3->x)
            {
                retstruct3->x = sum_val12;
                retstruct3->b1 = Z3->b1;
                retstruct3->a3 = CrossLayer::m_a3[i];
                retstruct3->zx = Z3->x;
                retstruct3->zy = Z3->y;
                retstruct3->zz = Z3->z;
            }
        }
        Z3 = Z3->next;
        //printf("\n");
    }
    //free(val_12);
    //if (retstruct3->x < 1)
        //cout << "when s1' = " << next_s1 << ", s2' = " << next_s2 << ", s3 = { " << curr_s3.x << ", " << curr_s3.y << "}, j = " << j << ", k = " << k << " Rin is " << Rin << ". optb1 is " << retstruct3->b1 << ", opta3 is " << retstruct3->a3 << ".\n";
    //cout << "value12 is " << retstruct3->x << ", optb1 is " << retstruct3->b1 << ", opta3 is " << retstruct3->a3 << ".\n";
    //printf("value12 is %f, optimal internal action is %d, optimal external action a3 is %d.\n", retstruct3->x, retstruct3->b1, retstruct3->a3);
    return retstruct3;
    
}

RETLY2 *stval1(float next_s1, float curr_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param)
{
    int i, j;
    float sum_val1 = 0.0;
    float prob2 = 1;        //CDMA
    RETLY2 *retstruct2 = new RETLY2;
    RETLY3 *val12 = new RETLY3;
    float thruwght, delwght, cstwght;
    
    thruwght = ThroughputWeight;
    delwght = DelayWeight;
    cstwght = CostWeight;
    
    retstruct2->x = 0.0;
    
    for(i = 0; i < 2; i++)   //iterate on a2
    {
        //sum_val1 = - cstwght * lamda2 * a2[i];
        sum_val1 = - Param.lamda2 * CrossLayer::m_a2[i];
        for(j = 0; j < 5; j++)  //for all next_s2, add value together
        {
            //printf("in stval1, i = %d, j = %d\n", i, j);
            //printf("next_s1 = %f in stval1.\n",next_s1);
            val12 = stval2(next_s1, CrossLayer::m_s2[j], Z3, curr_s3, NEXT_VAL, thruwght, delwght, cstwght, Param);
            //printf("val12 is %f for next_s2 = %f.\n",val12->x,next_s2[j]);
            /*if(i == 0)
            {
                if(s2[j] >= curr_s2)
                {
                    if(curr_s2 > 0.85)
                        prob2 = 0.1 / 2;
                    else
                        prob2 = 0.1 * 0.2 / (2 *(1 - curr_s2));
                }
                else
                    prob2 = 0.9 * 0.2 / (2 * (curr_s2 - 0.2));
                
            }
            else
            {
                if(s2[j] >= curr_s2)
                {
                    if(curr_s2 > 0.85)
                        prob2 = 0.9 / 2;
                    else
                        prob2 = 0.9 * 0.2 / (2 *(1 - curr_s2));
                }
                else
                    prob2 = 0.1 * 0.2 / (2 * (curr_s2 - 0.2));
            }*/
            prob2 = 0.2;
            
            sum_val1 += prob2 * val12->x;
        }
        retstruct2->a3 = val12->a3;
        retstruct2->b1 = val12->b1;
        retstruct2->zx = val12->zx;
        retstruct2->zy = val12->zy;
        retstruct2->zz = val12->zz;
        
        if(sum_val1 > retstruct2->x)
        {
            retstruct2->x = sum_val1;
            retstruct2->a2 = CrossLayer::m_a2[i];
        }
        
    }
    //cout << "value1 is " << retstruct2->x << ", optimal external action a2 is " << retstruct2->a2 << ".\n";
    //printf("value1 is %f, optimal external action a2 is %d.\n", retstruct2->x, retstruct2->a2);
    delete val12;
    return retstruct2;
}

RETLY1 *stvalfunc(float curr_s1, float curr_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param)
{
    int i, j, k, num;
    float sum_valcurr = 0.0;
    float tot_prob, ave_snr;
    float *prob1, *valx;
    float stateL;
    float *next_s1;
    RETLY1 *val_curr = new RETLY1;
    RETLY2 *val1 = new RETLY2;
    float thruwght, delwght, cstwght;
    float crossval, omval, ave_snrP, curr_s1P;
    
    thruwght = ThroughputWeight;
    delwght = DelayWeight;
    cstwght = CostWeight;
    
    /*std::cout << "possible external action a1s are ";
    for (i = 0; i < 10; i++)
        std::cout << CrossLayer::m_a1[i] << " ";
    std::cout << "\n";*/
    
    /*std::cout << "ThroughputWeight is " << thruwght << ", DelayWeight is " << delwght << ", CostWeight is " << cstwght << std::endl;
     std::cout << "PacketSize is " << Param.PktSize << ", PacketTime is " << Param.Tp << std::endl;*/
    
    val_curr->x = 0.0;
    //printf("curr_s1 is %f in stvalfunc.\n", curr_s1);
    if ((curr_s1 > 110) && (curr_s1) < 144)
    {
        num = 3;
        next_s1 = new float[3];
        next_s1[0] = curr_s1 - 2;
        next_s1[1] = curr_s1;
        next_s1[2] = curr_s1 + 2;
    }
    else if (curr_s1 < 110)
    {
        num = 2;
        next_s1 = new float[2];
        next_s1[0] = curr_s1;
        next_s1[1] = curr_s1 + 2;
        //printf("next_s1[0] = %f, next_s1[1] = %f\n", next_s1[0], next_s1[1]);
    }
    else
    {
        num = 2;
        next_s1 = new float[2];
        next_s1[0] = curr_s1 - 2;
        next_s1[1] = curr_s1;
    }

    for(i = 0; i < 10; i++)    //iterate on a1
    {
        //sum_valcurr = - cstwght * lamda1 * a1[i];
        sum_valcurr = - Param.lamda1 * CrossLayer::m_a1[i];
        //ave_snr = CrossLayer::m_a1[i] / 0.46;
        ave_snr = CrossLayer::m_a1[i] + 101;   //in dB
        ave_snrP = db2rt(ave_snr);
        curr_s1P = db2rt(curr_s1);
        
        //cout << "curr_s1 = " << curr_s1P << "\n";
        //cout << "curr_s1/ave_snr = " << curr_s1P/ave_snrP << "\n";
        //printf("average snr is %f when curr_s1 is %f.\n", ave_snrP, curr_s1P);
        //printf("average snr is %f when a1 is %f.\n", ave_snrP, CrossLayer::m_a1[i]);
        prob1 = new float[num];
        valx = new float[num];
        tot_prob = 0;
        
        for(j = 0; j < num; j++)
        {   //calculate p(s1'|s1, a1)
            if(next_s1[j] == curr_s1)
            {
                //if(next_s1[j] < 31)
                    //prob1[j] = 1 - (Param.Tp / omega(curr_s1, ave_snr)) * xing(0.4, ave_snr, Param.DopplerFreq);
                if(next_s1[j] > 144)
                {
                    crossval = xing(145, ave_snr, Param.DopplerFreq);
                    omval = omega(curr_s1, ave_snr);
                    prob1[j] = 1 - (Param.Tp / omega(curr_s1, ave_snr)) * crossval;
                    //cout << "value of xing is " << crossval << "\n";
                    //cout << "value of omega is " << omval << "\n";
                    //printf("prob of next_s1[%d] %f is %f, ", j, next_s1[j], prob1[j]);

                }
                else
                {
                    crossval = xing((next_s1[j] + 2), ave_snr, Param.DopplerFreq);
                    prob1[j] = 1 - (Param.Tp / omega(curr_s1, ave_snr)) * ( crossval - xing(next_s1[j], ave_snr, Param.DopplerFreq));
                }
            }
            else
            {
                stateL = max(curr_s1, next_s1[j]);
                crossval = xing(stateL, ave_snr, Param.DopplerFreq);
                prob1[j] = (Param.Tp / omega(curr_s1, ave_snr)) * crossval;
            }
            
            //cout << "value of xing is " << crossval << "\n";
            
            if(prob1[j] < 0)
                prob1[j] = 0;
            else
                tot_prob += prob1[j];
            
            //printf("in stvalfunc, i = %d, j = %d, num = %d\n", i, j, num);
            //printf("next_s1[j] = %f in stvalfunc.\n",next_s1[j]);
            //adding values for all next_s1
            val1 = stval1(next_s1[j], curr_s2, Z3, curr_s3, NEXT_VAL, thruwght, delwght, cstwght, Param);
            valx[j] = val1->x;
            
        }
        
        for (k = 0; k < num; k++)
        {
            prob1[k] = prob1[k] / tot_prob;
            sum_valcurr += prob1[k] * valx[k];
            //printf("prob of next_s1[%d] is %f, ", k, prob1[k]);
            
        }
        //printf("\n");
        

        val_curr->b1 = val1->b1;
        val_curr->a2 = val1->a2;
        val_curr->a3 = val1->a3;
        val_curr->zx = val1->zx;
        val_curr->zy = val1->zy;
        val_curr->zz = val1->zz;
        
        if(sum_valcurr > val_curr->x)
        {
            val_curr->x = sum_valcurr;
            val_curr->a1 = CrossLayer::m_a1[i];
        }
        
    }
    //cout << "current value is " << val_curr->x << ", optimal external action a1 is " << val_curr->a1 << ".\n";
    //printf("current value is %f, optimal external action a1 is %d.\n", val_curr->x, val_curr->a1);
    delete val1;
    delete next_s1;
    return val_curr;
}


/*RETLY1 *stvalfunc(float PHY, float MAC, APPST APP, QOS *Z3, float *NEXT_VAL)
{
    RETLY1 *valcurr = (RETLY1 *)malloc(sizeof(RETLY1));
    
    valcurr = stvalcurr(PHY, s1, s2, Z3, APP, NEXT_VAL);
    
    return valcurr;
}*/

