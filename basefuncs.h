//
//  basefuncs.h
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#ifndef xlayer_basefuncs_h
#define xlayer_basefuncs_h

int fact(int x);
QOS *add_to_list(float x, float y, float z, int b1, QOS **list);
void deallocmem(QOS *list);

#endif
