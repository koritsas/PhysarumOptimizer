package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

/**
 * Created by ilias on 5/10/2016.
 */

public abstract class PhysarumPolycephalumSP extends AbstractPhysarumPolycephalum {
    public PhysarumPolycephalumSP(Graph graph, Node sourceNode,Node sinkNode,double Io,double γ,int numberOfIterations){
        this.graph=graph;
        this.sourceNode=sourceNode;
        this.sinkNode=sinkNode;
        this.numberOfIterations=numberOfIterations;
        this.γ=γ;
        this.Io=Io;
    }

    @Override
    public double calculateTubeConductivity(Edge e) {
        double D = conductivityMap.get(e);

        double Q = fluxMap.get(e);
        double L=getEdgeCost(e);

        double p1 =pressureMap.get(e.getNodeA());
        double p2 =pressureMap.get(e.getNodeB());

        double ps=pressureMap.get(sourceNode);
        double pe =pressureMap.get(sinkNode);

         double newD =(0.5)*((Q*(p1-p2))/(L*(ps-pe))+D);

        //double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));
        //double fQ=Math.abs(Q);

       // double newD = fQ - 0.1 * D;



        return newD;
    }
}
