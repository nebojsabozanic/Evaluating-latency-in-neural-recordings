        /*
 * SPIKE-distance (Kreuz) doi:
 *
 * @Nebojsa, @Mario April 1 2015
 *
 * 2 do:
 *      add overall simlarity - Done
 *      add 1 inline for nearest neigbour
 *      free memory within the loops - Done (i hope so)
 *      rename variables
 *      add a fool proof *tEnd > *tBeg, and some more !!!
 */

#include "mex.h"

double abso(double a) {
    return (a > 0) ? a : -a;
}

template <class T> const T& min (const T& a, const T& b) {
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

template <class T> const T& max (const T& a, const T& b) {
  return !(b>a)?a:b;
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
#define spike_distance_out plhs[0]
    
#define spikes_in prhs[0]
#define tBeg_in prhs[1]
#define tEnd_in prhs[2]
#define rBeg_in prhs[3]
#define rEnd_in prhs[4]
    
    int nTrains, nPairs, trainCount[2], spikeCount[2], nSpikes[2], cc, lcc, llcc, nnindx[2], flag[2], pac = 0, debug = 0;
    double *spike_distance_profile, *tBeg, *tEnd, t, t_prev, *spkDist, *spikes[2], sLeft, sRight, isiLeft[2], isiRight[2], S[2], nearestPrevious[2], nearestFuture[2], *rBeg, *rEnd;
    const mxArray *spikesPr[2];
    
    tBeg = mxGetPr(tBeg_in);
    tEnd = mxGetPr(tEnd_in);
    rBeg = mxGetPr(rBeg_in);
    rEnd = mxGetPr(rEnd_in);
    
    /* protections for rEnd > max(spikes) > min(spikes) > rBeg, rEnd > tEnd > tBeg > rBeg*/
    nTrains = mxGetN(spikes_in);
    nPairs = nTrains * (nTrains - 1)/2;
    
    spike_distance_out = mxCreateNumericArray(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetM(spike_distance_out, 1);
    mxSetN(spike_distance_out, nPairs);
    mxSetData(spike_distance_out, mxMalloc(nPairs*sizeof(double)));
    spkDist = mxGetPr(spike_distance_out);
    
    if (*tEnd < *tBeg)
        mexPrintf("ERROR; Time frame ill defined");
    else {
        
        for(lcc = 0; lcc < nPairs; ++lcc)
            spkDist[lcc] = 0;
        
        for(trainCount[0] = 0; trainCount[0] < nTrains-1; ++trainCount[0]) {
            spikesPr[0] = mxGetCell(spikes_in, trainCount[0]);
            spikes[0] = mxGetPr(spikesPr[0]);
            nSpikes[0] = mxGetN(spikesPr[0]);
            
            if (nSpikes[0] > 0) {
                
                for(trainCount[1] = trainCount[0] + 1; trainCount[1] < nTrains;  ++trainCount[1]) {
                    spikesPr[1] = mxGetCell(spikes_in, trainCount[1]);
                    spikes[1] = mxGetPr(spikesPr[1]);
                    nSpikes[1] = mxGetN(spikesPr[1]);
                    
                    if (nSpikes[1] > 0) {
                        
                        t = *tBeg;
                        spikeCount[0] = 0; spikeCount[1] = 0;
                        nnindx[0] = 0; nnindx[1] = 0;
                        flag[0] = 0; flag[1] = 0;
                        
                        /* ######################################## Before time tBeg ############################################# */
                        for(lcc = 0; lcc < 2; ++lcc) {
                            
                            if (nSpikes[lcc] > 1) {
                                while (t >= spikes[lcc][spikeCount[lcc]+1] && spikeCount[lcc] < nSpikes[lcc]-1)
                                    ++spikeCount[lcc];
                            }
                            if (t >= spikes[lcc][spikeCount[lcc]])
                                isiRight[lcc] = ((spikeCount[lcc] < nSpikes[lcc]-1) ? spikes[lcc][spikeCount[lcc]+1] : rEnd[lcc])- spikes[lcc][spikeCount[lcc]];
                            else {
                                isiRight[lcc] = max(spikes[lcc][min(1, nSpikes[lcc-1])] - spikes[lcc][0], spikes[lcc][0] - rBeg[lcc]);
								
                                while (abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][min(nnindx[lcc]+1, nSpikes[1-lcc]-1)]) < abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][nnindx[lcc]]))
                                    ++nnindx[lcc];
                                
                                nearestPrevious[lcc] = abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][nnindx[lcc]]);
                                
                                /*auxilary? a bit vague theory...*/
                                if (nnindx[lcc] == nSpikes[1-lcc] -1)
                                    if (nearestPrevious[lcc] > abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]]))
                                        nearestPrevious[lcc] = abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]]);
                                if (nnindx[lcc] == 0)
                                    if (nearestPrevious[lcc] > abso(spikes[lcc][spikeCount[lcc]] - rBeg[1-lcc]))
                                        nearestPrevious[lcc] = abso(spikes[lcc][spikeCount[lcc]] - rBeg[1-lcc]);
                                
                                if (t < spikes[lcc][spikeCount[lcc]] || spikeCount[lcc] == nSpikes[lcc] - 1)
                                    nearestFuture[lcc] = nearestPrevious[lcc];
                                else {
                                    while (abso(spikes[lcc][spikeCount[lcc]+1] - spikes[1-lcc][min(nnindx[lcc]+1, nSpikes[1-lcc]-1)]) < abso(spikes[lcc][spikeCount[lcc]+1] - spikes[1-lcc][nnindx[lcc]]))
                                        ++nnindx[lcc];
                                    nearestFuture[lcc] = abso(spikes[lcc][spikeCount[lcc]+1] - spikes[1-lcc][nnindx[lcc]]);
                                    
                                    if (nnindx[lcc] == nSpikes[1-lcc] -1)
                                        if (nearestFuture[lcc] > abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]+1]))
                                            nearestFuture[lcc] = abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]+1]);
                                    if (nnindx[lcc] == 0)
                                        if (nearestFuture[lcc] > abso(spikes[lcc][spikeCount[lcc]+1] - rBeg[1-lcc]))
                                            nearestFuture[lcc] = abso(spikes[lcc][spikeCount[lcc]+1] - rBeg[1-lcc]);
                                }
                            }
                        }
                        
                        for(cc = 0; cc < 2; ++cc) {
							if (spikes[cc][spikeCount[cc]] <= t)
									S[cc] = (nearestPrevious[cc] * (((spikeCount[cc] < nSpikes[cc]-1) ? spikes[cc][spikeCount[cc]+1] : rEnd[cc]) - t) + nearestFuture[cc]*(t - spikes[cc][spikeCount[cc]]))/(((spikeCount[cc] < nSpikes[cc]-1) ? spikes[cc][spikeCount[cc]+1] : rEnd[cc]) - spikes[cc][spikeCount[cc]]); /*max*/
							else
                                S[cc] = (nearestPrevious[cc] * (spikes[cc][spikeCount[cc]] - t) + nearestFuture[cc]*(t - rBeg[cc]))/(spikes[cc][spikeCount[cc]] - rBeg[cc]);
                        
                            if (spikeCount[cc] == 0 && spikes[cc][0] > t) 
                                if (spikes[cc][0] - rBeg[cc]< ((nSpikes[cc] > 1 ? spikes[cc][1] : rEnd[cc]) - spikes[cc][0]))
									if (nSpikes[cc] > 1)
										S[cc] *= (spikes[cc][0] - rBeg[cc])/(spikes[cc][1] - spikes[cc][0]);
                            
                            if (spikes[cc][spikeCount[cc]] <= t) {
                                if (spikeCount[cc] < nSpikes[cc] - 1)
                                    ++spikeCount[cc];
                                else
                                    flag[cc] = 1;
                            }
                        }
                        
                        sRight = 2 * (S[0] * isiRight[1] + S[1] * isiRight[0]) / ((isiRight[0] + isiRight[1])*(isiRight[0] + isiRight[1])); /*fix t- for conor*/
                        /*mexPrintf("%f\n", sRight);*/
                        
                        /* ######################################## Between time tBeg and tEnd ############################################# */
                        while ((!flag[0] || !flag[1]) && t < *tEnd) {
                          
                            t_prev = t;
                            
                            if (spikes[0][spikeCount[0]] == spikes[1][spikeCount[1]]) {
                                t = spikes[0][spikeCount[0]];
                                for(lcc = 0; lcc < 2; ++lcc)
                                    if (spikeCount[lcc] == nSpikes[lcc]-1)
                                        flag[lcc] = 1;
                                }
                            else {
                                if (flag[1])
                                    cc = 0; 
                                else {
                                    if(flag[0])
                                        cc = 1;
                                    else {
                                        if (spikes[0][spikeCount[0]] < spikes[1][spikeCount[1]]) 
                                            cc = 0;
                                        else
                                            cc = 1;
                                    } 
                                }
                                if (spikeCount[cc] == nSpikes[cc]-1)
                                    flag[cc] = 1;
                                t = spikes[cc][spikeCount[cc]];
                            }
                            
                            if (t < *tEnd) {
                                
                                if (spikes[0][spikeCount[0]] == spikes[1][spikeCount[1]]) {
                                    
                                    for(lcc = 0; lcc < 2; ++lcc) {
                                        if (spikeCount[lcc] < nSpikes[lcc]-1)
                                            isiRight[lcc] = spikes[lcc][spikeCount[lcc]+1] - spikes[lcc][spikeCount[lcc]];
                                        else {
											if (nSpikes[lcc] > 1)
												isiRight[lcc] = max(rEnd[lcc] - spikes[lcc][nSpikes[lcc]-1], spikes[lcc][nSpikes[lcc]-1] - spikes[lcc][nSpikes[lcc]-2]);
											else
												isiRight[lcc] = rEnd[lcc] - spikes[lcc][nSpikes[lcc]-1];

                                            nearestPrevious[lcc] = 0;
                                            
                                            if (spikeCount[lcc] < nSpikes[lcc]-1)
                                                ++spikeCount[lcc];
                                            
                                            while (abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][min(nnindx[lcc]+1, nSpikes[1-lcc]-1)]) < abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][nnindx[lcc]]))
                                                ++nnindx[lcc];
                                            
                                            /*auxilary? a bit vague theory...*/
                                            nearestFuture[lcc] = abso(spikes[lcc][spikeCount[lcc]] - spikes[1-lcc][nnindx[lcc]]); /*if same ;)*/
                                            if (nnindx[lcc] == nSpikes[1-lcc] -1)
                                                if (nearestFuture[lcc] > abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]]))
                                                    nearestFuture[lcc] = abso(rEnd[1-lcc] - spikes[lcc][spikeCount[lcc]]);
                                            if (nnindx[lcc] == 0)
                                                if (nearestFuture[lcc] > abso(spikes[lcc][spikeCount[lcc]] - rBeg[1-lcc]))
                                                    nearestFuture[lcc] = abso(spikes[lcc][spikeCount[lcc]] - rBeg[1-lcc]);
                                        }
                                    }
                                    
                                    S[0] = 0;
                                    S[1] = 0;
                                    spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] = spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] + sRight/2 * (t - t_prev);
                                
                                    sRight = 0;
                                    
                                    /*mexPrintf("%f\n", sRight);
                                    mexPrintf("%f\n", sRight);*/
                                }
                                else {
                              
                                    t = spikes[cc][spikeCount[cc]];
                                    
                                    isiLeft[cc] = isiRight[cc]; /* within the loop*/
                                    
                                    if (spikeCount[cc] < nSpikes[cc]-1)
                                        isiRight[cc] = spikes[cc][spikeCount[cc]+1] - spikes[cc][spikeCount[cc]];
                                    else {
										if (nSpikes[cc] > 1)
											isiRight[cc] = max(rEnd[cc] - spikes[cc][nSpikes[cc]-1], spikes[cc][nSpikes[cc]-1] - spikes[cc][nSpikes[cc]-2]);
										else
											isiRight[cc] = rEnd[cc] - spikes[cc][nSpikes[cc]-1];

                                        isiLeft[1-cc] = isiRight[1-cc];
                                        
                                        nearestPrevious[cc] = nearestFuture[cc];
                                      
                                        if (spikeCount[cc] < nSpikes[cc]-1)
                                            ++spikeCount[cc];

                                        while (abso(spikes[cc][spikeCount[cc]] - spikes[1-cc][min(nnindx[cc]+1, nSpikes[1-cc]-1)]) < abso(spikes[cc][spikeCount[cc]] - spikes[1-cc][nnindx[cc]]))
                                            ++nnindx[cc];
                                        nearestFuture[cc] = abso(spikes[cc][spikeCount[cc]] - spikes[1-cc][nnindx[cc]]);
                                        
                                        /*auxilary? a bit vague theory...*/
                                        if (nnindx[cc] == nSpikes[1-cc] -1)
                                            if (nearestFuture[cc] > abso(rEnd[1-cc] - spikes[cc][spikeCount[cc]]))
                                                nearestFuture[cc] = abso(rEnd[1-cc] - spikes[cc][spikeCount[cc]]);
                                        if (nnindx[cc] == 0)
                                            if (nearestFuture[cc] > abso(spikes[cc][spikeCount[cc]] - rBeg[1-cc]))
                                                nearestFuture[cc] = abso(spikes[cc][spikeCount[cc]] - rBeg[1-cc]);
                                        
                                        if (spikeCount[cc] > 1)
                                            S[cc] = nearestPrevious[cc];
                                            
                                        if (spikes[1-cc][spikeCount[1-cc]] <= t) {
                                            if (spikeCount[1-cc] < nSpikes[1-cc]-1 && spikeCount[1-cc] > 0) 
                                                 S[1-cc] = (nearestPrevious[1-cc] * (spikes[1-cc][spikeCount[1-cc]+1] - t) + nearestFuture[1-cc]*(t - (spikes[1-cc][spikeCount[1-cc]])))/isiLeft[1-cc];
                                            }
                                        else {
                                            if (nSpikes[1-cc] > 1)
                                                S[1-cc] = (nearestPrevious[1-cc] * (spikes[1-cc][spikeCount[1-cc]] - t) + nearestFuture[1-cc]*(t - ((spikeCount[1-cc] > 0) ? spikes[1-cc][spikeCount[1-cc]-1] : rBeg[1-cc])))/isiLeft[1-cc];
                                            
                                            sLeft = 2 * (S[0] * isiLeft[1] + S[1] * isiLeft[0]) / ((isiLeft[0] + isiLeft[1])*(isiLeft[0] + isiLeft[1]));
                                            /*mexPrintf("%f\n", sLeft);*/
                                            S[cc] = nearestPrevious[cc];
                                            
                                            spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] = spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] + (sLeft + sRight)/2 * (t - t_prev);
                                            
                                            if (flag[cc])
                                                if (rEnd[cc] - spikes[cc][nSpikes[cc]-1] < (spikes[cc][nSpikes[cc]-1] - (nSpikes[cc] > 1 ? spikes[cc][nSpikes[cc]-2] : rBeg[cc])))
													if (nSpikes[cc] > 1)
														S[cc] *= (rEnd[cc] - spikes[cc][nSpikes[cc]-1])/(spikes[cc][nSpikes[cc]-1] - spikes[cc][nSpikes[cc]-2]);
                                             
                                            sRight = 2 * (S[0] * isiRight[1] + S[1] * isiRight[0]) / ((isiRight[0] + isiRight[1])*(isiRight[0] + isiRight[1]));
                                            /*mexPrintf("%f\n", sRight);*/
                                        }
                                    }
                                }
                            }
                        }
                        
                        /* ######################################## After time tEnd ############################################# */
                        
                        if (t < *tEnd)
                            t_prev = t;
                        
                        t = *tEnd;
                        
                        isiLeft[0] = isiRight[0];
                        isiLeft[1] = isiRight[1];
                        
                        for(cc = 0; cc < 2; ++cc) {
                            if (spikes[cc][spikeCount[cc]] < t) {
                                if (spikeCount[cc] < nSpikes[cc]-1)
                                    S[cc] = (nearestPrevious[cc] * (spikes[cc][spikeCount[cc]+1] - t) + nearestFuture[cc]*(t - (spikes[cc][spikeCount[cc]])))/isiLeft[cc];
                            }
                            else {
                                if ((nSpikes[cc] > 1) && (spikeCount[cc] == nSpikes[cc]-1))
                                    S[cc] = (nearestPrevious[cc] * (spikes[cc][spikeCount[cc]] - t) + nearestFuture[cc]*(t - ((spikeCount[1-cc] > 0) ? spikes[1-cc][spikeCount[1-cc]-1] : rBeg[1-cc])))/isiLeft[cc];
                        
                                sLeft = 2 * (S[0] * isiLeft[1] + S[1] * isiLeft[0]) / ((isiLeft[0] + isiLeft[1])*(isiLeft[0] + isiLeft[1])); /*fix t- for conor*/
                                /*mexPrintf("%f\n", sLeft);*/
                                spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] = spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] + (sLeft + sRight)/2 * (t - t_prev);
                                
                                spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1] = spkDist[trainCount[0]*(trainCount[0]+1)/2 + trainCount[1]-1]/(*tEnd - *tBeg);
                            }
                        }       
                    }
                }
            }
        }
    }
    return;
}

/* free the memory, delete  mxDestroyArray  mxFree*/