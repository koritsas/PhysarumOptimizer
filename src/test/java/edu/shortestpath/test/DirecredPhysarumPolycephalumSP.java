package edu.shortestpath.test;

import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicDirectedNode;

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
        double c=0;
        if(node1.getInEdge(node2)!=null){
            c=-super.calculateCoefficient(node1,node2);
        }else{
            c=super.calculateCoefficient(node1,node2);
        }

        return c;

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
