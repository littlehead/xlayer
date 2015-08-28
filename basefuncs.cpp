//
//  otherfuncs.c
//  xlayer
//
//  Created by Yuping Dong on 9/20/13.
//  Copyright (c) 2013 Yuping Dong. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <new>
#include "xlayer-modulel.hpp"
#include "types.hpp"
#include "basefuncs.hpp"

using namespace std;

int fact(int x)
{
    while (x != 0)
        return x * fact(x - 1);
    return 1;
}


// adding data to the end of list, returning list head
QOS *add_to_list(float x, float y, float z, int b1, QOS **list)
{
    QOS *listtemp;
    QOS *listtemp1;
    
    listtemp1 = new QOS;
    //printf("Creating list failed.\n");
    //return NULL;
    
    listtemp1->x = x;
    listtemp1->y = y;
    listtemp1->z = z;
    listtemp1->b1 = b1;
    listtemp1->next = NULL;
    
    
    if(!(*list))
    {
        //printf("adding first element...\n");
    
        *list = listtemp1;
    }
    else
    {
        //printf("adding next element...\n");
        listtemp = *list;
        while(listtemp->next != NULL)
            listtemp = listtemp->next;
        listtemp->next = listtemp1;
    }
    
    return *list;
}

void deallocmem(QOS *list)
{
    QOS *temp = list;
    QOS *node = NULL;
    while (temp != NULL)        //wipes off memory block
    {
        node = temp;
        temp = temp->next;
        delete node;
    }
    //list = NULL;                //clears pointer  why can't I combine these two into one function?
}                                 //list is now local variable, so this only changes the local one. In main, Z1 is
                                  //still not changed.
