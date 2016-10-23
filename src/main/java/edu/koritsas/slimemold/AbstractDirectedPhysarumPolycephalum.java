package edu.koritsas.slimemold;

import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.DirectedEdge;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicDirectedNode;

import java.util.List;

/**
 * Created by ilias on 21/10/2016.
 */
public abstract class AbstractDirectedPhysarumPolycephalum extends AbstractPhysarumPolycephalum {
    public AbstractDirectedPhysarumPolycephalum(DirectedGraph graph, Node sourceNode,Node sinkNode,double Io,double γ,int numberOfIterations){
        super(graph,sourceNode,sinkNode,Io,γ,numberOfIterations);
        this.graph=graph;

    }
    @Override
    public double calculateTubeConductivity(Edge e) {
        DirectedEdge edge= (DirectedEdge) e;
        double D = conductivityMap.get(edge);


        double Q = fluxMap.get(e);


        double L = getEdgeCost(e);



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
    @Override
    protected double calculateCoefficient(Node n1, Node n2) {
        BasicDirectedNode node1 = (BasicDirectedNode) n1;
        BasicDirectedNode node2 = (BasicDirectedNode) n2;
        List<Edge> ndf =node1.getInEdges();

        double c=0;
        Edge edge =node1.getEdge(node2);
        if(node1.getOutEdges().contains(edge)){

            c=super.calculateCoefficient(node1,node2);
        }else if(node1.getInEdges().contains(edge)){
            c=-FastMath.abs(super.calculateCoefficient(node1,node2));
        }else{
            c=0;
        }


        return c;

    }

    @Override
    protected double calculateTubeFlux(Edge e) {
        DirectedEdge edge = (DirectedEdge) e;
        double p1 =pressureMap.get(edge.getInNode());
        double p2 =pressureMap.get(edge.getOutNode());

        double D =conductivityMap.get(edge);
        double w = getEdgeCost(e);

        double Q=(D/w)*(p1-p2);
        if(Q<0){
            Q=0;
        }


        return Q;

    }
}
