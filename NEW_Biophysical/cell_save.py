import pickle
import saveClass as sc
import libcell as lb
import numpy as np
import struct
import os

# def save_Ldend(Ldends, bfname):
#     # create a binary file
#     bfname='Dend_length.bin'
#     binfile = open(bfname, 'wb')
#     # and write out two integers with the row and column dimension
#     header = struct.pack('2I', Ldends.shape[0], Ldends.shape[1])
#     binfile.write(header)
#     # then loop over columns and write each
#     for i in range(Ldends.shape[1]):
#         ddata = struct.pack('%id' % Ldends.shape[0], *Ldends[:,i])
#         binfile.write(ddata)
#     binfile.close()

def save_ave_replay(aveData, nIter, nStart, bfname):
    vd = np.zeros((nIter, 4, nStart))

    for i_trial in range(nIter):
        vv = aveData[i_trial]
        for i_dendrite in range(4):
            vvv = vv[i_dendrite]
            mv = np.reshape(vvv, (nStart, 1501))
            vd[i_trial, i_dendrite, :] = np.mean(mv[:,550:1000], 1)

    mvd = np.mean(vd, 0)

    # print (bfname)

    # create a binary file
    binfile = open(bfname, 'wb')
    # and write out two integers with the row and column dimension
    header = struct.pack('2I', mvd.shape[0], mvd.shape[1])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mvd.shape[1]):
        ddata = struct.pack('%id' % mvd.shape[0], *mvd[:,i])
        binfile.write(ddata)
    binfile.close()

def save_ave_place(aveData, nIter, bfname):
    vd = np.zeros((nIter, 4, 20))

    for i_trial in range(nIter):
        vv = aveData[i_trial]
        for i_dendrite in range(4):
            vvv = vv[i_dendrite]
            mv = np.reshape(vvv[0:50000], (20, 2500))
            vd[i_trial, i_dendrite, :] = np.mean(mv, 1)

    mvd = np.mean(vd, 0)

    print (bfname)

    # create a binary file
    binfile = open(bfname, 'wb')
    # and write out two integers with the row and column dimension
    header = struct.pack('2I', mvd.shape[0], mvd.shape[1])
    binfile.write(header)
    # then loop over columns and write each
    for i in range(mvd.shape[1]):
        ddata = struct.pack('%id' % mvd.shape[0], *mvd[:,i])
        binfile.write(ddata)
    binfile.close()


def save_sim(data, out_binary=True, out_vdend=True, out_pickle=False, outdir='data', dt_save=1.0):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    modelData = sc.emptyObject()
    lb.props(modelData)

    if (data.stimType=='DStim'):
        filename = 'T' + str(data.TSTOP) + '_dend' + str(data.iclampLoc[2]) + '_N' + str(len(data.iRange)) + '_I' + str(data.iRange[0]) + '_dI' + str(data.iRange[1]-data.iRange[0]) 
    elif (data.stimType=='SStim'):
        filename = 'T' + str(data.TSTOP) + '_soma_N' + str(len(data.iRange)) + '_I' + str(data.iRange[0]) + '_dI' + str(data.iRange[1]-data.iRange[0]) 

    else :
        filename = 'T' + str(data.TSTOP) + '_Ne' + str(data.Ensyn)+'_gA'+str(round(data.Agmax,2)) + '_tauA' + str(data.Atau2)
        if (data.NMDA):
            filename = filename + '_gN'+str(round(data.Ngmax,2))
        if (data.GABA):
            filename = filename  + '_Ni'+str(data.Insyn) + '_gG'+str(round(data.Igmax, 2))
        if (data.GABA_B):
            filename = filename  + '_gB'+str(round(data.Bgmax, 2))

        if (data.modulateNa):
            filename = filename  + '_noDendNa'

        if (data.stimType == 'nIter'):
            filename = filename + '_tInt' + str(data.tInterval) + 'ms_' + data.locBias + '_' + data.direction
        
        if (data.stimType == 'place'):
            filename = filename + "_Er" + str(data.Erate) + '_Ir'+str(data.Irate) + '_' + data.placeType + '_rep' + str(data.nIter)
            filename = filename + '_stimseed' + str(data.stimseed)

    if (data.randomW == True):
        filename = filename + '_randW'

    if out_pickle:
        dataList = [data, modelData]
        fname = './'+outdir+'/'+filename+'.pkl'
        f = open(fname, 'wb')
        pickle.dump(dataList, f)
        f.close()


    if out_binary:
        #---------------------------------------------
        # WRITE the response in a binary file to read it with R
        mat = np.array(data.vdata)
        L = mat.shape[1]
        #dt_ratio = int(round(dt_save / data.dt))
        dt_ratio = int(dt_save / data.dt)
        print(dt_ratio)
        mat = mat[:,0:L:dt_ratio]

        #bfname = './'+outdir+'/vdata_'+filename+'.bin'
        bfname = './'+outdir+'/vdata_'+filename+'.npy'
        print (bfname)
        np.save(bfname, mat)

        """
        # create a binary file
        binfile = open(bfname, 'wb')
        # and write out two integers with the row and column dimension
        header = struct.pack('2I', mat.shape[0], mat.shape[1])
        binfile.write(header)
        # then loop over columns and write each
        for i in range(mat.shape[1]):
            ddata = struct.pack('%id' % mat.shape[0], *mat[:,i])
            binfile.write(ddata)
        binfile.close()
        """

        if out_vdend:
            # WRITE the dendritic response
            nRep = len(data.vDdata)
            mat = np.array(data.vDdata[0])
            for i in range(1, nRep):
                mat = np.hstack((mat, data.vDdata[i]))

            L = mat.shape[1]
            dt_ratio = int(round(dt_save / data.dt))
            mat = mat[:,0:L:dt_ratio]
            
            #bfname = './'+outdir+'/vDdata_'+filename+'.bin'
            bfname = './'+outdir+'/vDdata_'+filename+'.npy'
            np.save(bfname, mat)

            """
            # create a binary file
            binfile = open(bfname, 'wb')
            # and write out two integers with the row and column dimension
            header = struct.pack('2I', mat.shape[0], mat.shape[1])
            binfile.write(header)
            # then loop over columns and write each
            for i in range(mat.shape[1]):
                ddata = struct.pack('%id' % mat.shape[0], *mat[:,i])
                binfile.write(ddata)
            binfile.close()
            """

        # ---------------------------------------------
        # WRITE the location of the synapses        
        if (data.GABA) :
            Ilocs = np.array(data.Ilocs) 
            Ilocs[:,1] = 1 + Ilocs[:,1] # code that these are inhibitory synapses
            Elocs = np.array(data.Elocs)
            Locs = np.row_stack((Elocs, Ilocs))
        else :
            Locs = np.array(data.Elocs)

        #bfname = './'+outdir+'/synlocs_'+filename+'.bin'
        bfname = './'+outdir+'/synlocs_'+filename+'.npy'
        print (bfname)
        np.save(bfname, Locs)

        """
        # create a binary file
        binfile = open(bfname, 'wb')
        # and write out two integers with the row and column dimension
        header = struct.pack('2I', Locs.shape[0], Locs.shape[1])
        binfile.write(header)
        # then loop over columns and write each
        for i in range(Locs.shape[1]):
            ddata = struct.pack('%id' % Locs.shape[0], *Locs[:,i])
            binfile.write(ddata)
        binfile.close()
        """

        #---------------------------------------------
        # Write the input spike train
        if (len(data.stim)>0):
            stim = data.stim
            #bfname = './'+outdir+'/stim_'+filename+'.bin'
            bfname = './'+outdir+'/stim_'+filename+'.npy'

            np.save(bfname, stim)
            
            """
            # create a binary file
            binfile = open(bfname, 'wb')
            # and write out two integers with the row and column dimension
            header = struct.pack('2I', stim.shape[0], stim.shape[1])
            binfile.write(header)
            # then loop over columns and write each
            for i in range(stim.shape[1]):
                ddata = struct.pack('%id' % stim.shape[0], *stim[:,i])
                binfile.write(ddata)
            binfile.close()
            """
