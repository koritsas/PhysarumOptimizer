package edu.shortestpath.test;

import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.structure.DirectedEdge;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicDirectedNode;

import java.util.List;

/**
 * Created by ilias on 12/10/2016.
 */
public abstract class DirecredPhysarumPolycephalumSP extends AbstractPhysarumPolycephalum {

    public DirecredPhysarumPolycephalumSP(DirectedGraph graph, Node source, Node sink, double Io, double γ, int numberOfIterations){
        this.graph=graph;
        this.sourceNode=source;
        this.sinkNode=sink;
        this.Io=Io;
        this.γ=γ;
        this.numberOfIterations =numberOfIterations;
    }
/*
    @Override
    protected RealVector createConstantsVector() {
        List<Node> allNodes =getAllNodes(graph);
        double[] constants = new double[allNodes.size()];
        for (int i = 0; i < constants.length; i++) {

            Node n1 = allNodes.get(i);
            if (n1.equals(sourceNode)) {

                constants[i] = Io + calculateCoefficient(n1,sinkNode) * pressureMap.get(sinkNode);

            } else if (n1.equals(sinkNode)) {
                constants[i] = -Io - calculateSelfCoefficient(n1) *  pressureMap.get(sinkNode);
            } else {

                constants[i] = calculateCoefficient(n1,sinkNode) *  pressureMap.get(sinkNode);


            }
        }


        return new ArrayRealVector(constants);
    }
*/
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

    /*
        @Override
        public double calculateSelfCoefficient(Node n) {
            BasicDirectedNode node = (BasicDirectedNode) n;


            ToDoubleFunction<Edge> function = new ToDoubleFunction<Edge>() {
                @Override
                public double applyAsDouble(Edge edge) {

                    return calculateCoefficient(edge.getNodeA(),edge.getNodeB());
                }
            };

            List<Edge> inEdges =node.getInEdges();
            List<Edge> outEdges =node.getOutEdges();
          double  c=inEdges.stream().collect(Collectors.summingDouble(function)).doubleValue()-outEdges.stream().collect(Collectors.summingDouble(function)).doubleValue();

            return c;
        }
    */
    @Override
    public double calculateTubeConductivity(Edge e) {
        DirectedEdge edge= (DirectedEdge) e;
        double D = conductivityMap.get(edge);


        double Q = fluxMap.get(e);
        double newD;

            double L = getEdgeCost(e);

            double p1 = pressureMap.get(edge.getInNode());
            double p2 = pressureMap.get(edge.getOutNode());

            double ps = pressureMap.get(sourceNode);
            double pe = pressureMap.get(sinkNode);



            newD = (0.5) * ((Q * (p1 - p2)) / (L * (ps - pe)) + D);



        //double fQ = Math.pow(Math.abs(Q), γ) / (1 + Math.pow(Math.abs(Q), γ));
       // double fQ=Math.abs(Q);

        // double newD = fQ - 0.4 * D;



        return newD;
    }
}
