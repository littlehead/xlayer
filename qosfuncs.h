//
//  qosfuncs.h
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#ifndef xlayer_qosfuncs_h
#define xlayer_qosfuncs_h

QOS *optqosfrontier(QOS *QOS_lower, float State, int ActionB, float ThroughputWeight, float DelayWeight, float CostWeight);
QOS *calcphyqos(float PhyState, int ModLevel, float PktTime, int PktSize);
void calcqosphy(float PhyState, int ModLevel, float PktTime, int PktSize, QOS **QOSPhy);

#endif
