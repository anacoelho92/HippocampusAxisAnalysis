#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 09:33:45 2019

@author: anacoelho

Script to improve hippocampus centerline (previously computed), 
define orthogonal planes and compute volumes of sections

Usage: 
  centerline_analysis <working_dir> <subjects_file> <roi_name> <LoR> <nplanes>
  
  working_dir: directory with hippocampus 3D stl file, centerline points and 
               control points
  subjects_file: file with subjects codes
  roi_name: hippo (same name in other files)
  LoR: lh - left; rh - right
  nplanes: number of planes to define (usually 20)
    
"""

import numpy as np
import vtk
import sys
import os


# function to close rendering window
def close_window(iren):
    render_window = iren.GetRenderWindow()
    render_window.Finalize()
    iren.TerminateApp()

# function to find the closest point to a specific set of coordinates
def get_closest_point_fromCoords(point,coords):
    min_id = -1
    min_dist = 9999999
    for i in range(len(coords)):
        p = coords[i]
        dist = np.sqrt(pow(p[0]-point[0],2)+pow(p[1]-point[1],2)+pow(p[2]-point[2],2))
        if dist < min_dist:
            min_dist = dist
            min_id = i

    return min_id

# function to compute areas of the orthogonal planes
def compute_areas(nplanes,centerline_coords,plane_centers, surface, renderer,display):
    

    planes_col = vtk.vtkPlaneCollection() #collection to save all planes
    plane_centers = np.array(plane_centers) #coordinates of the center of each plane
    
    # create planes (using center and normal vector)
    for i in range(nplanes):

        # center
        coords1 = plane_centers[i,:] 
        
        # get the point in the centerline that is closest to this plane center 
        j = get_closest_point_fromCoords(coords1,centerline_coords)

        # define two points (before and after the closest point) in the centerline to compute the normal vector of the plane
        coords2 = centerline_coords[j-1,:]   
        coords3 = centerline_coords[j+1,:]
        
        # compute normal and origin of the plane
        norm = (coords3[0]-coords2[0],coords3[1]-coords2[1],coords3[2]-coords2[2])
        orig = coords1
        
        # set normal and origin of the plane and add it to the collection
        plane=vtk.vtkPlane()
        plane.SetOrigin(orig)
        plane.SetNormal(norm)
        planes_col.AddItem(plane)
        
    areas=[]
    centerOfMass = []
    vectors = []
    
    # cut the surface with the defined planes and compute areas
    for i in range(planes_col.GetNumberOfItems()):
        plane = planes_col.GetItem(i)
        vectors.append(plane.GetNormal()) # save normal vector of all planes
                
        # cut the surface
        cutter = vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInputConnection(surface.GetOutputPort())
        cutter.Update()
            
        # get the contour line cutting the surface
        stripper = vtk.vtkStripper()
        stripper.SetInputConnection(cutter.GetOutputPort())
        stripper.Update()
        circle = stripper.GetOutput()
                
        # transform line to polydata
        boundaryPoly = vtk.vtkPolyData()
        boundaryPoly.SetPoints(circle.GetPoints())
        boundaryPoly.SetPolys(circle.GetLines())
            
        # triangulate the polygon to create cut surface
        triangles = vtk.vtkTriangleFilter()
        triangles.SetInputData(boundaryPoly)
        triangles.Update()
                
        # get area of cut surface
        mass=vtk.vtkMassProperties()
        mass.SetInputData(triangles.GetOutput())
        area_cut = mass.GetSurfaceArea()
        areas.append(area_cut)
        
        # compute center of mass of cut surface
        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInputData(triangles.GetOutput())
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()
        centerOfMass.append(centerOfMassFilter.GetCenter())
        
        # create actors to display in the scene the contour line and the cut surface
        circleMapper = vtk.vtkPolyDataMapper()
        circleMapper.SetInputData(circle)
        circleActor = vtk.vtkActor()
        circleActor.SetMapper(circleMapper)
        circleActor.GetProperty().SetColor(1,0,0)
                
        boundaryMapper = vtk.vtkPolyDataMapper()
        boundaryMapper.SetInputData(triangles.GetOutput())
        boundaryActor = vtk.vtkActor()
        boundaryActor.SetMapper(boundaryMapper)
        boundaryActor.GetProperty().SetColor(1,1,0)
        
        # add the actors to the scene 
        if display==1:
            renderer.AddActor(circleActor)
            renderer.AddActor(boundaryActor)
        
    return areas,centerOfMass,planes_col,vectors


# function to compute volumes of the sections
def compute_volumes(planes_col, surface, renderer, display):

    volumes = [] # array to save volumes of all sections
    
    # extremities are special cases
    # 1st plane
    cut_planes_col1 = vtk.vtkPlaneCollection() # plane used to cut section
        
    p1 = planes_col.GetItem(0)
    plane1=vtk.vtkPlane()
    plane1_norm = p1.GetNormal()
    plane1.SetOrigin(p1.GetOrigin())
    plane1.SetNormal((-plane1_norm[0],-plane1_norm[1],-plane1_norm[2]))
    cut_planes_col1.AddItem(plane1)
     
    # clip surface with plane
    clipper1 = vtk.vtkClipPolyData()
    clipper1.SetInputData(surface.GetOutput())
    clipper1.SetClipFunction(plane1)
    clipper1.Update()
           
    # create a cap
    boundaryEdges1 = vtk.vtkFeatureEdges()
    boundaryEdges1.SetInputConnection(clipper1.GetOutputPort())
    boundaryEdges1.BoundaryEdgesOn()
    boundaryEdges1.FeatureEdgesOff()
    boundaryEdges1.NonManifoldEdgesOff()
    
    boundaryClean1 = vtk.vtkCleanPolyData()
    boundaryClean1.SetInputConnection(boundaryEdges1.GetOutputPort())
    
    boundaryStrips1 = vtk.vtkStripper()
    boundaryStrips1.SetInputConnection(boundaryClean1.GetOutputPort())
    boundaryStrips1.Update()
    
    boundaryPoly1 = vtk.vtkPolyData()
    boundaryPoly1.SetPoints(boundaryStrips1.GetOutput().GetPoints())
    boundaryPoly1.SetPolys(boundaryStrips1.GetOutput().GetLines())
    
    boundaryTriangles1 = vtk.vtkTriangleFilter()
    boundaryTriangles1.SetInputData(boundaryPoly1)
    boundaryTriangles1.Update()
    
    # add clipped surface and the cap
    appendFilter1 = vtk.vtkAppendPolyData()
    appendFilter1.AddInputData(clipper1.GetOutput())
    appendFilter1.AddInputData(boundaryTriangles1.GetOutput())
    appendFilter1.Update()
    
    cleanFilter1 = vtk.vtkCleanPolyData()
    cleanFilter1.SetInputConnection(appendFilter1.GetOutputPort())
    cleanFilter1.Update()
  
    normals1 = vtk.vtkPolyDataNormals()
    normals1.SetInputData(cleanFilter1.GetOutput())
    normals1.ConsistencyOn()
    normals1.SplittingOff()
    normals1.Update()
  
    # check for open edges
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.SetInputConnection(normals1.GetOutputPort())
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.Update()
#    print("Number of open edges = ",featureEdges.GetOutput().GetNumberOfCells())
    
    # compute volume
    mass1 = vtk.vtkMassProperties()
    mass1.SetInputConnection(normals1.GetOutputPort())
    mass1.Update()
    volumes.append(mass1.GetVolume())
            
    # create actor to display section in the scene
    clipMapper1 = vtk.vtkDataSetMapper()
    clipMapper1.SetInputConnection(normals1.GetOutputPort())
    
    clipActor1 = vtk.vtkActor()
    clipActor1.SetMapper(clipMapper1)
    clipActor1.GetProperty().SetInterpolationToFlat()
    color_rgb = list(np.random.choice(range(256), size=3))    
    color = [x / 256 for x in color_rgb]
    clipActor1.GetProperty().SetColor(color)
    
    
    if display==1:
        renderer.AddActor(clipActor1)
    
    
    # intermediate planes
    for i in range(1,planes_col.GetNumberOfItems()):

        # collection with planes used to cut section
        cut_planes_col = vtk.vtkPlaneCollection()
        
        # first plane
        p0 = planes_col.GetItem(i-1)
        plane0=vtk.vtkPlane()
        plane0.SetOrigin(p0.GetOrigin())
        plane0.SetNormal(p0.GetNormal())
        cut_planes_col.AddItem(plane0)
        
        # second plane
        p1 = planes_col.GetItem(i)
        plane1=vtk.vtkPlane()
        plane1_norm = p1.GetNormal()
        plane1.SetOrigin(p1.GetOrigin())
        plane1.SetNormal((-plane1_norm[0],-plane1_norm[1],-plane1_norm[2]))
        cut_planes_col.AddItem(plane1)
        
        # clip surface with the two planes
        clipper0 = vtk.vtkClipPolyData()
        clipper0.SetInputData(surface.GetOutput())
        clipper0.SetClipFunction(plane0)
        clipper0.Update()
        
        clipper = vtk.vtkClipPolyData()
        clipper.SetInputConnection(clipper0.GetOutputPort())
        clipper.SetClipFunction(plane1)
        clipper.Update()
        
        # create a cap  
        # cut the surface
        cutter0 = vtk.vtkCutter()
        cutter0.SetCutFunction(plane0)
        cutter0.SetInputConnection(surface.GetOutputPort())
        cutter0.Update()
            
        # get the line cutting the surface
        stripper0 = vtk.vtkStripper()
        stripper0.SetInputConnection(cutter0.GetOutputPort())
        stripper0.Update()
        circle0 = stripper0.GetOutput()
                
        # transform line to polydata
        boundaryPoly0 = vtk.vtkPolyData()
        boundaryPoly0.SetPoints(circle0.GetPoints())
        boundaryPoly0.SetPolys(circle0.GetLines())
            
        # triangulate the polygon to create cut surface
        boundaryTriangles0 = vtk.vtkTriangleFilter()
        boundaryTriangles0.SetInputData(boundaryPoly0)
        boundaryTriangles0.Update()
        
        # cut the surface
        cutter1 = vtk.vtkCutter()
        cutter1.SetCutFunction(plane1)
        cutter1.SetInputConnection(surface.GetOutputPort())
        cutter1.Update()
            
        # get the line cutting the surface
        stripper1 = vtk.vtkStripper()
        stripper1.SetInputConnection(cutter1.GetOutputPort())
        stripper1.Update()
        circle1 = stripper1.GetOutput()
                
        # transform line to polydata
        boundaryPoly1 = vtk.vtkPolyData()
        boundaryPoly1.SetPoints(circle1.GetPoints())
        boundaryPoly1.SetPolys(circle1.GetLines())
            
        # triangulate the polygon to create cut surface
        boundaryTriangles1 = vtk.vtkTriangleFilter()
        boundaryTriangles1.SetInputData(boundaryPoly1)
        boundaryTriangles1.Update()
        
        
        # add clipped surface and the cap
        appendFilter = vtk.vtkAppendPolyData()
        appendFilter.AddInputData(clipper.GetOutput())
        appendFilter.AddInputData(boundaryTriangles0.GetOutput())
        appendFilter.AddInputData(boundaryTriangles1.GetOutput())
        appendFilter.Update()
        
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
        cleanFilter.Update()
    
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(cleanFilter.GetOutput())
        normals.ConsistencyOn()
        normals.SplittingOff()
        normals.Update()
        
        boundaryEdges = vtk.vtkFeatureEdges()
        boundaryEdges.SetInputConnection(normals.GetOutputPort())
        boundaryEdges.BoundaryEdgesOn()
        boundaryEdges.FeatureEdgesOff()
        boundaryEdges.NonManifoldEdgesOff()
        
        boundaryClean = vtk.vtkCleanPolyData()
        boundaryClean.SetInputConnection(boundaryEdges.GetOutputPort())
        
        boundaryStrips = vtk.vtkStripper()
        boundaryStrips.SetInputConnection(boundaryClean.GetOutputPort())
        boundaryStrips.Update()
        
        boundaryPoly = vtk.vtkPolyData()
        boundaryPoly.SetPoints(boundaryStrips.GetOutput().GetPoints())
        boundaryPoly.SetPolys(boundaryStrips.GetOutput().GetLines())
        
        boundaryTriangles = vtk.vtkTriangleFilter()
        boundaryTriangles.SetInputData(boundaryPoly)
        boundaryTriangles.Update()
        
        
        # add clipped surface and the cap
        appendFilter1 = vtk.vtkAppendPolyData()
        appendFilter1.AddInputData(normals.GetOutput())
        appendFilter1.AddInputData(boundaryTriangles.GetOutput())
        appendFilter1.Update()
        
        cleanFilter1 = vtk.vtkCleanPolyData()
        cleanFilter1.SetInputConnection(appendFilter1.GetOutputPort())
        cleanFilter1.Update()
    
        normals1 = vtk.vtkPolyDataNormals()
        normals1.SetInputData(cleanFilter1.GetOutput())
        normals1.ConsistencyOn()
        normals1.SplittingOff()
        normals1.Update()
        
        # check for open edges
        featureEdges = vtk.vtkFeatureEdges()
        featureEdges.SetInputConnection(normals1.GetOutputPort())
        featureEdges.FeatureEdgesOff()
        featureEdges.BoundaryEdgesOn()
        featureEdges.NonManifoldEdgesOn()
        featureEdges.Update()
#        print("Number of open edges = ",featureEdges.GetOutput().GetNumberOfCells())
    
        # compute volume
        mass = vtk.vtkMassProperties()
        mass.SetInputConnection(normals1.GetOutputPort())
        mass.Update()
        volumes.append(mass.GetVolume())
        
        # create actor to display section in the scene
        clipMapper = vtk.vtkDataSetMapper()
        clipMapper.SetInputConnection(normals1.GetOutputPort())
            
        clipActor = vtk.vtkActor()
        clipActor.SetMapper(clipMapper)
        clipActor.GetProperty().SetInterpolationToFlat()
        color_rgb = list(np.random.choice(range(256), size=3))
        color = [x / 256 for x in color_rgb]
        clipActor.GetProperty().SetColor(color)
        
            
        if display==1:
            renderer.AddActor(clipActor)
            
    # last plane
    cut_planes_col2 = vtk.vtkPlaneCollection() # plane used to cut section
        
    p2 = planes_col.GetItem(planes_col.GetNumberOfItems()-1)
    plane2=vtk.vtkPlane()
    plane2.SetOrigin(p2.GetOrigin())
    plane2.SetNormal(p2.GetNormal())
    cut_planes_col2.AddItem(plane2)
        
    # clip surface with plane
    clipper2 = vtk.vtkClipPolyData()
    clipper2.SetInputData(surface.GetOutput())
    clipper2.SetClipFunction(plane2)
    clipper2.Update()
            
    # create a cap
    boundaryEdges2 = vtk.vtkFeatureEdges()
    boundaryEdges2.SetInputConnection(clipper2.GetOutputPort())
    boundaryEdges2.BoundaryEdgesOn()
    boundaryEdges2.FeatureEdgesOff()
    boundaryEdges2.NonManifoldEdgesOff()
    
    boundaryClean2 = vtk.vtkCleanPolyData()
    boundaryClean2.SetInputConnection(boundaryEdges2.GetOutputPort())
    
    boundaryStrips2 = vtk.vtkStripper()
    boundaryStrips2.SetInputConnection(boundaryClean2.GetOutputPort())
    boundaryStrips2.Update()
    
    boundaryPoly2 = vtk.vtkPolyData()
    boundaryPoly2.SetPoints(boundaryStrips2.GetOutput().GetPoints())
    boundaryPoly2.SetPolys(boundaryStrips2.GetOutput().GetLines())
    
    boundaryTriangles2 = vtk.vtkTriangleFilter()
    boundaryTriangles2.SetInputData(boundaryPoly2)
    boundaryTriangles2.Update()
    
    # add clipped surface and the cap
    appendFilter2 = vtk.vtkAppendPolyData()
    appendFilter2.AddInputData(clipper2.GetOutput())
    appendFilter2.AddInputData(boundaryTriangles2.GetOutput())
    appendFilter2.Update()
    
    cleanFilter2 = vtk.vtkCleanPolyData()
    cleanFilter2.SetInputConnection(appendFilter2.GetOutputPort())
    cleanFilter2.Update()
    
    normals2 = vtk.vtkPolyDataNormals()
    normals2.SetInputData(cleanFilter2.GetOutput())
    normals2.ConsistencyOn()
    normals2.SplittingOff()
    normals2.Update()
    
    # check for open edges
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.SetInputConnection(normals2.GetOutputPort())
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.Update()
#    print("Number of open edges = ",featureEdges.GetOutput().GetNumberOfCells())
    
    # compute volume
    mass2 = vtk.vtkMassProperties()
    mass2.SetInputConnection(normals2.GetOutputPort())
    mass2.Update()
    volumes.append(mass2.GetVolume())
    
    # create actor to display section in the scene
    clipMapper2 = vtk.vtkDataSetMapper()
    clipMapper2.SetInputConnection(normals2.GetOutputPort())
    
    clipActor2 = vtk.vtkActor()
    clipActor2.SetMapper(clipMapper2)
    clipActor2.GetProperty().SetInterpolationToFlat()
    color_rgb = list(np.random.choice(range(256), size=3))   
    color = [x / 256 for x in color_rgb]
    clipActor2.GetProperty().SetColor(color)
    
    if display==1:
        renderer.AddActor(clipActor2)
    
    return volumes
 
# function to create line between two points
def addLine(renderer, p1, p2, color=[0.0, 0.0, 1.0]):
    
    # create line object
    line = vtk.vtkLineSource()
    line.SetPoint1(p1)
    line.SetPoint2(p2)

    # create actor and add to scene
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(line.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    renderer.AddActor(actor)
   
    
# function to create a point
def addPoint(renderer, p, radius=0.25, color=[0.0, 0.0, 0.0]):
    
    # create point object
    point = vtk.vtkSphereSource()
    point.SetCenter(p)
    point.SetRadius(radius)
    point.SetPhiResolution(100)
    point.SetThetaResolution(100)

    # create actor and add to scene
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(point.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    renderer.AddActor(actor)
  
# function to improve bezier line
def interpolate_bezier(start_pt,end_pt,centers,nplanes,renderer,display):
    
    # planes centers where the line should pass
    centers = np.array(centers)
    
    # intermediate point between first and last of the line
    midpt=[0,0,0]
    midpt[0] = (end_pt[0] + start_pt[0])/2
    midpt[1] = (end_pt[1] + start_pt[1])/2
    midpt[2] = (end_pt[2] + start_pt[2])/2
    
    # get the plane center that is closest to the intermediate point
    id_mid = get_closest_point_fromCoords(midpt,centers)
    id_mid = int(nplanes / 2)
    
    # define control points for the bezier line using the 3 points (first, intermediate and last)
    ctrl_pt=[0,0,0]
    ctrl_pt[0] = 2 * centers[id_mid][0] - 0.5 * start_pt[0] - 0.5 * end_pt[0]
    ctrl_pt[1] = 2 * centers[id_mid][1] - 0.5 * start_pt[1] - 0.5 * end_pt[1]
    ctrl_pt[2] = 2 * centers[id_mid][2] - 0.5 * start_pt[2] - 0.5 * end_pt[2]
    
    # number of points to define the bezier line
    nsamples = 1000
    t = np.linspace(0,1,1000)
    bezier_pts = np.zeros((nsamples,3))
    
    # compute bezier line point coordinates using bezier line equations
    for i in range(0,nsamples):
        B=[0,0,0]
        
        B[0] = (1-t[i])*(1-t[i])*start_pt[0] + 2*t[i]*(1-t[i])*ctrl_pt[0] + t[i]*t[i]*end_pt[0]
        B[1] = (1-t[i])*(1-t[i])*start_pt[1] + 2*t[i]*(1-t[i])*ctrl_pt[1] + t[i]*t[i]*end_pt[1]
        B[2] = (1-t[i])*(1-t[i])*start_pt[2] + 2*t[i]*(1-t[i])*ctrl_pt[2] + t[i]*t[i]*end_pt[2]
    
        bezier_pts[i,:] = B


    # create object with points
    new_centerline_pts = vtk.vtkPoints()   
    
    for i in range(len(bezier_pts)):
        new_centerline_pts.InsertNextPoint(bezier_pts[i][0],bezier_pts[i][1],bezier_pts[i][2])
       
    # create line object
    new_polyLine = vtk.vtkPolyLine()
    new_polyLine.GetPointIds().SetNumberOfIds(new_centerline_pts.GetNumberOfPoints())
                
    for i in range(new_centerline_pts.GetNumberOfPoints()):
        new_polyLine.GetPointIds().SetId(i,i)
                
    new_cells = vtk.vtkCellArray()
    new_cells.InsertNextCell(new_polyLine)
        
    new_centerline_pd = vtk.vtkPolyData()
    new_centerline_pd.SetPoints(new_centerline_pts)
    new_centerline_pd.SetLines(new_cells)
    
    # create actor for the line
    new_centerlineMapper = vtk.vtkPolyDataMapper()
    new_centerlineMapper.SetInputData(new_centerline_pd)
        
    new_centerlineActor = vtk.vtkActor()
    new_centerlineActor.SetMapper(new_centerlineMapper)
    new_centerlineActor.GetProperty().SetColor(1,1,1)
    
    # get equidistant points from line that will be the new plane centers
    new_pts=np.zeros((nplanes,3))
    percent = 1/(nplanes+1)
    for r in range(1,(nplanes+1)):
        
        d = (percent * r)
        x = (1-d) * (1-d) * start_pt[0] + 2 * (1-d) * d * ctrl_pt[0] + d * d * end_pt[0]
        y = (1-d) * (1-d) * start_pt[1] + 2 * (1-d) * d * ctrl_pt[1] + d * d * end_pt[1]
        z = (1-d) * (1-d) * start_pt[2] + 2 * (1-d) * d * ctrl_pt[2] + d * d * end_pt[2]
        
        new_pts[r-1,0] = x
        new_pts[r-1,1] = y
        new_pts[r-1,2] = z
        
    
    # add equidistant points, first and last line points and line to the scene
    if display==1:
        for i in range(len(new_pts)):
            addPoint(renderer,new_pts[i,:],color=[0.0, 0.0, 1.0])
        addPoint(renderer,start_pt,color=[0,0,1])
        addPoint(renderer,end_pt,color=[0,0,1])
        renderer.AddActor(new_centerlineActor)
        
         
    return bezier_pts, new_pts
 

# main function that corrects the line, computes planes and volumes of sections
def centerline_analysis(roi_file, centerline_file,ctrl_pts_file, nplanes, LoR, out_dir):
    
    # load hippocampus STL file    
    reader_stl = vtk.vtkSTLReader()
    reader_stl.SetFileName(roi_file)
    reader_stl.Update()

    # smooth surface
    smoothFilter = vtk.vtkWindowedSincPolyDataFilter()
    smoothFilter.SetInputConnection(reader_stl.GetOutputPort())
    smoothFilter.SetNumberOfIterations(25)
    smoothFilter.SetPassBand(0.001)
    smoothFilter.NonManifoldSmoothingOn()
    smoothFilter.NormalizeCoordinatesOn()
    smoothFilter.Update()
    
    # check for open edges
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.SetInputConnection(smoothFilter.GetOutputPort())
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.Update()
#    print("Number of open edges = ",featureEdges.GetOutput().GetNumberOfCells())
    
    # get total volume 
    mass = vtk.vtkMassProperties()
    mass.SetInputConnection(smoothFilter.GetOutputPort())
    mass.Update()
    total_vol = mass.GetVolume()
    
    # create actor to add surface to the scene
    surfMapper = vtk.vtkPolyDataMapper()
    surfMapper.SetInputConnection(smoothFilter.GetOutputPort())
    surfMapper.ScalarVisibilityOff()
        
    surfActor = vtk.vtkActor()
    surfActor.SetMapper(surfMapper)
    surfActor.GetProperty().SetOpacity(0.2)
    surfActor.GetProperty().SetColor(217/255,217/255,217/255)
        
    # rendering things ...
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(.5,.5,.5)
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    
    # create line from coordinates (calculated previously)
    centerline_coords = np.loadtxt(centerline_file)
    
    centerline_coords=sorted(centerline_coords, key=lambda tup: (-tup[2]))
    centerline_coords=np.asarray(centerline_coords)
    
    # load control points (calculated previously)
    ctrl_pts = np.loadtxt(ctrl_pts_file)
    ctrl_pts=sorted(ctrl_pts, key=lambda tup: (-tup[2]))
    ctrl_pts=np.asarray(ctrl_pts)

    # create object for line points
    centerline_pts = vtk.vtkPoints()   
    
    for i in range(len(centerline_coords)):
        centerline_pts.InsertNextPoint(centerline_coords[i][0],centerline_coords[i][1],centerline_coords[i][2])
        
    centerline_pd = vtk.vtkPolyData()
    centerline_pd.SetPoints(centerline_pts)
    
    
    # get first and last point of the line that intersect with the surface
    # compute intersection between surface and line points
    selectEnclosedPoints = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPoints.SetInputData(centerline_pd)
    selectEnclosedPoints.SetSurfaceData(smoothFilter.GetOutput())
    selectEnclosedPoints.Update()
    
    first=[]
    last=[]
    # iterate over first half of the line and get the first point of intersection
    for i in range(int(len(centerline_coords)/2)):
        if selectEnclosedPoints.IsInside(i)==1:
            first = i
            break
        
    # iterate over last half of the line and get the last point of intersection
    for i in range(len(centerline_coords)-1,int(len(centerline_coords)/2),-1):
        if selectEnclosedPoints.IsInside(i)==1:
            last = i
            break

    # redefine centerline coordinates
    centerline_coords=centerline_coords[first:last,:]

    # update points object with redefined points
    centerline_pts = vtk.vtkPoints()   
    
    for i in range(len(centerline_coords)):
        centerline_pts.InsertNextPoint(centerline_coords[i][0],centerline_coords[i][1],centerline_coords[i][2])
       
    # create object for the line
    polyLine = vtk.vtkPolyLine()
    polyLine.GetPointIds().SetNumberOfIds(centerline_pts.GetNumberOfPoints())
                
    for i in range(centerline_pts.GetNumberOfPoints()):
        polyLine.GetPointIds().SetId(i,i)
                
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyLine)
        
    centerline_pd = vtk.vtkPolyData()
    centerline_pd.SetPoints(centerline_pts)
    centerline_pd.SetLines(cells)
    
    # create actor to add line to the scene
    centerlineMapper = vtk.vtkPolyDataMapper()
    centerlineMapper.SetInputData(centerline_pd)
        
    centerlineActor = vtk.vtkActor()
    centerlineActor.SetMapper(centerlineMapper)
    centerlineActor.GetProperty().SetColor(0,1,0)

        
    # compute t value for end points (relative position in the line)
    a = ctrl_pts[0][0] - (2 * ctrl_pts[1][0]) + ctrl_pts[2][0]
    b = -(2 * ctrl_pts[0][0]) + (2 * ctrl_pts[1][0])
    c0 = ctrl_pts[0][0] - centerline_coords[0][0]
    c1 = ctrl_pts[0][0] - centerline_coords[len(centerline_coords)-1][0]
    coeff0 = [a, b, c0]
    results0 = np.roots(coeff0)
    coeff1 = [a, b, c1]
    results1 = np.roots(coeff1)
    
    t0 = min(results0)
    t1 = max(results1)
   
    
    # compute points for the center of the planes
    new_pts=np.zeros((nplanes,3))
    
    percent = (t1-t0)/(nplanes+1) # spacing between points
   
    for r in range(1,(nplanes+1)):

        t = (percent * r) + t0 # relative position in the line
        # compute coordinates using line equations
        x = (1-t) * (1-t) * ctrl_pts[0][0] + 2 * (1-t) * t * ctrl_pts[1][0] + t * t * ctrl_pts[2][0]
        y = (1-t) * (1-t) * ctrl_pts[0][1] + 2 * (1-t) * t * ctrl_pts[1][1] + t * t * ctrl_pts[2][1]
        z = (1-t) * (1-t) * ctrl_pts[0][2] + 2 * (1-t) * t * ctrl_pts[1][2] + t * t * ctrl_pts[2][2]
        
        new_pts[r-1,0] = x
        new_pts[r-1,1] = y
        new_pts[r-1,2] = z
          
    
    # create planes and compute area of intersection of each plane
    areas,centerOfMass, planes_col = compute_areas(nplanes,centerline_coords,new_pts, smoothFilter, renderer,0)
    
    # update line to pass through the center of mass of all planes
    new_centerline_coords, new_centers = interpolate_bezier(centerline_coords[0],centerline_coords[len(centerline_coords)-1],centerOfMass,nplanes,renderer,1)

    # create new planes with the updated line
    areas,centerOfMass, planes_col = compute_areas(nplanes,new_centerline_coords,new_centers, smoothFilter, renderer,1)

    # save centers and normals of all planes to create nifti images for tractography
    np.savetxt(os.path.join(out_dir, LoR + '_centers.txt'), centerOfMass)
    np.savetxt(os.path.join(out_dir, LoR + '_normals.txt'), planes_vectors)
    
    # add first and last line points to the scene with different colors
    addPoint(renderer,new_centerline_coords[0,:], color=[1,0,0])
    addPoint(renderer,new_centerline_coords[len(new_centerline_coords)-1,:], color=[0,1,0])
 
    # clip surface with planes and compute volumes
    volumes=[]
    volumes = compute_volumes(planes_col, smoothFilter, renderer,0)
    
    # add actors to scene and render it (line and surface)
    renderer.AddActor(surfActor)
    renderer.AddActor(centerlineActor)
    renderWindow.Render()
    renderWindowInteractor.Start()

    
    close_window(renderWindowInteractor)
    del renderWindow, renderWindowInteractor
    
    return (areas, volumes, total_vol)



# run the main function
# check if we have all arguments
if len(sys.argv) < 6:
    print ('usage: centerline_analysis <working_dir> <subjects_file> <roi_name> <LoR> <nplanes>')
else:
    # process inputs
    working_dir = str(sys.argv[1]) # working directory
    subjects_filepath = str(sys.argv[2]) # file with subjects ID
    roi_name = str(sys.argv[3]) # ROI name
    LoR = str(sys.argv[4]) # Left or Right
    nplanes = int(sys.argv[5]) # number of planes

    # iterate over subjects
    with open(subjects_filepath, 'r') as subjects:
        mylist = subjects.read().splitlines()
        #print(len(mylist))
        areas = np.zeros((nplanes,len(mylist)))
        volumes = np.zeros((nplanes+1,len(mylist)))
        total_volumes = np.zeros((len(mylist),1))
        
        for i in range(len(mylist)):
            
            # get subject name from file
            subject = mylist[i]
            print (subject)
            
            # hippocampus STL file
            roi_sub_path = os.path.join(working_dir,subject, subject+'_'+LoR+'_'+roi_name+'.stl')
            # centerline coordinates file
            centerline_sub_path = os.path.join(working_dir, subject,subject+'_'+LoR+'_'+roi_name+'_bezier_pts.txt')
            # control points file
            ctrl_pts_path = os.path.join(working_dir, subject,subject+'_'+LoR+'_'+roi_name+'_ctrl_pts.txt')
            # run centerline analysis
            out_dir = os.path.join(working_dir, subject) # directory to save planes center and normal vector (to create nifti images for tractography)
            subj_areas, subj_vol, subj_total_vol = centerline_analysis(roi_sub_path, centerline_sub_path, ctrl_pts_path, nplanes, LoR, out_dir)
            
            # save areas and volumes of all subjects
            areas[:,i] = subj_areas
            volumes[:,i] = subj_vol
            total_volumes[i,0] = subj_total_vol

    # save results
    np.savetxt(os.path.join(working_dir, 'volumes_' + LoR + '.txt'), volumes)
    np.savetxt(os.path.join(working_dir, 'areas_' + LoR + '.txt'), areas)
    np.savetxt(os.path.join(working_dir, 'totalVolumes_' + LoR + '.txt'), total_volumes) 
