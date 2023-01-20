import networkx as nx
import vtk
import numpy as np

from config import simInputData


def save_VTK(sid, G, apertures, fracture_lens, lens, flow, pressure, cb, name): 
    """
    This function is written using VTK module
    Input:
        1. Network graph G
        2. isolated_nodes from "Remove_isolated_nodes" method to avoid the data error by fillinf zeros there
        3. Suffix (string) to save name of files

    """

    np={} 

    pos = list(zip(nx.get_node_attributes(G, 'x').values(), nx.get_node_attributes(G, 'y').values(), nx.get_node_attributes(G, 'z').values()))


        # Generate the polyline for the spline. 
    points = vtk.vtkPoints() 
    edgeData = vtk.vtkPolyData() 

        # Edges 

    lines = vtk.vtkCellArray()          
        

    point_data = vtk.vtkDoubleArray()
    point_data.SetNumberOfComponents(2)
    point_data.SetComponentName(0, 'p')
    point_data.SetComponentName(1, 'cb')
    for n in G.nodes():
        (x,y,z) = pos[n]
        points.InsertPoint(n,x,y,z)
        point_data.InsertNextTuple([pressure[n], cb[n]])
    
    #Filling zeros at deleted nodes
    # try:
    #     for i in self.isolated_nodes:
    #         points.InsertPoint(int(i),0,0,0)
    # except:
    #     pass
    
    cell_data_d = vtk.vtkDoubleArray()
    cell_data_d.SetNumberOfComponents(1)
    cell_data_d.SetName('b')
    
    cell_data_l = vtk.vtkDoubleArray()
    cell_data_l.SetNumberOfComponents(1)
    cell_data_l.SetName('l')

    cell_data_fl = vtk.vtkDoubleArray()
    cell_data_fl.SetNumberOfComponents(1)
    cell_data_fl.SetName('fl')

    cell_data_q = vtk.vtkDoubleArray()
    cell_data_q.SetNumberOfComponents(1)
    cell_data_q.SetName('q')
    
    cell_data_cb_in = vtk.vtkDoubleArray()
    cell_data_cb_in.SetNumberOfComponents(1)
    cell_data_cb_in.SetName('cb in')
   
    cell_data_cb_out = vtk.vtkDoubleArray()
    cell_data_cb_out.SetNumberOfComponents(1)
    cell_data_cb_out.SetName('cb_out')

    tmp_u = []; tmp_v = [];
    for i, e in enumerate(G.edges()):
        u=e[0] 
        v=e[1]
        lines.InsertNextCell(2)  
        lines.InsertCellPoint(u) 
        lines.InsertCellPoint(v)
        cell_data_d.InsertNextTuple([apertures[i] * sid.b0])
        cell_data_l.InsertNextTuple([lens[i] * sid.l0])
        cell_data_fl.InsertNextTuple([fracture_lens[i] * sid.l0])
        cell_data_q.InsertNextTuple([flow[i]])
        cell_data_cb_in.InsertNextTuple([cb[u]])
        cell_data_cb_out.InsertNextTuple([cb[v]])


    edgeData.GetCellData().AddArray(cell_data_d)
    edgeData.GetCellData().AddArray(cell_data_l)
    edgeData.GetCellData().AddArray(cell_data_fl)
    edgeData.GetCellData().AddArray(cell_data_q)
    edgeData.GetCellData().AddArray(cell_data_cb_in)
    edgeData.GetCellData().AddArray(cell_data_cb_out)


    edgeData.SetPoints(points) 
    edgeData.SetLines(lines)

    writer = vtk.vtkXMLPolyDataWriter();
    writer.SetFileName(sid.dirname + '/' + name);

    writer.SetInputData(edgeData)
    writer.Write()

def load_VTK(file_name):
    
    """
    
        This function reads the Network model output VTK files
        ==========================
        "Without Pore-Merging"
        ==========================
        and returns a Dictionary containing points and scalars in the file
        
        INPUT:
            VTK file saved during simulations
           
        OUTPUT:
            A Dictionary containing cordinates of network points and all scalar fields
               
    """
    #Dictionary where the data of VTK file
    scalars = {}
    
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(file_name)
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    points = np.array(data.GetPoints().GetData())
    pointdata = data.GetPointData()
    celldata = data.GetCellData()

    #Extracting Cell bounds or pore nodes
    CellArray = data.GetCells()
    pores = CellArray.GetData()
    no_of_pores = CellArray.GetNumberOfCells()
    
    pores_array = np.array(pores)
    list_of_pore_connection = []

    for i in range(no_of_pores):
        list_of_pore_connection.append([pores_array[j] for j in range(i*3, i*3+3)])
    scalars['Cell_Nodes'] = list_of_pore_connection;
    
    #No of properties saved in VTK files
    no_of_fields = reader.GetNumberOfScalarsInFile() 

    scalars['Points'] = points;
    for i in range(no_of_fields): 
        name_of_field = reader.GetScalarsNameInFile(i)
        scalars[name_of_field] = np.array(celldata.GetArray(name_of_field))
                
    return scalars


def build_VTK(sid:simInputData):

    scalars = load_VTK(sid.vtk_name)

    n = sid.n
    nkw = n ** 2
    l = sid.l

    in_nodes = list(range(n))
    out_nodes = list(range(nkw - n, nkw))
    G_edges = []
    diams = []
    lens = []
    boundary_edges = []
    n0 = 0
    d0 = 0
    for i,e in enumerate(scalars['Cell_Nodes']):
        d = scalars['d'][i]
        if d > 0:
            d0 += d
            n0 += 1

    d0 = d0 / n0
    print (d0)
    
    sid.d0 = d0
 
    for i,e in enumerate(scalars['Cell_Nodes']):
        n1 = e[1]
        n2 = e[2]
        d = scalars['d'][i] / sid.d0
        q = scalars['q'][i]
        
        l = 1 #0.730316
        if d != 0:
            if not ((n1 < n and n2 >= nkw - n) or (n2 < n and n1 >= nkw - n)):
                G_edges.append((n1, n2))
                diams.append(d)
                lens.append(l)

            if (n1 % n == n-1 and n2 % n == 0) or (n2 % n == n-1 and n1 % n == 0):
                boundary_edges.append((n1, n2))
    
    nodes = []
    for node in scalars['Points']:
        nodes.append((node[0], node[1]))

    G = nx.Graph()
    G.add_nodes_from(list(range(nkw)))

    G.add_edges_from(G_edges)

    for node in G.nodes:
        G.nodes[node]["pos"] = nodes[node]

    for i, edge in enumerate(G.edges()):
        G[edge[0]][edge[1]]['d'] = diams[i]
        G[edge[0]][edge[1]]['l'] = lens[i]

    return G, in_nodes, out_nodes, boundary_edges
    