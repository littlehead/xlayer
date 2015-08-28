//
//  iterfuncs.h
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#ifndef xlayer_iterfuncs_h
#define xlayer_iterfuncs_h

/*#ifndef max
#define max (a, b) ((a > b) ? (a) : (b))
#endif*/


RETLY3 *stval2(float next_s1, float next_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param);
RETLY2 *stval1(float next_s1, float curr_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param);
RETLY1 *stvalfunc(float curr_s1, float curr_s2, QOS *Z3, APPST curr_s3, float NEXT_VAL[][5][25], float ThroughputWeight, float DelayWeight, float CostWeight, ValueIterationParameters Param);
//RETLY1 *stvalfunc(float PHY, float MAC, APPST APP, QOS *Z3, float *NEXT_VAL);

#endif
