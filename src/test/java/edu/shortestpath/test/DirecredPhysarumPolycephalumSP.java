package edu.shortestpath.test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.DirectedNode;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Node;

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

    @Override
    protected RealMatrix createCoefficientsMatrix() {
        List<Node> allNodes=getAllNodes(graph);
        List<Node> allButSink =getAllButSink(graph);

        double[][] matrix = new double[allNodes.size()][allButSink.size()];

        for (int i = 0; i < allNodes.size(); i++) {
            DirectedNode n1 = (DirectedNode) allNodes.get(i);
            for (int j = 0; j < allButSink.size(); j++) {
                DirectedNode n2 = (DirectedNode) allButSink.get(j);

                if(n1.equals(n2)){
                    matrix[i][j]=calculateSelfCoefficient(n1);
                }else{
                    if(n1.getOutEdge(n2)!=null){
                        matrix[i][j]=-calculateCoefficient(n1.getOutEdge(n2));
                    }else if (n1.getInEdge(n2)!=null){
                        matrix[i][j]=calculateCoefficient(n1.getInEdge(n2));
                    }else{
                        matrix[i][j]=0;
                    }

                }


            }
        }


        return new Array2DRowRealMatrix(matrix);
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
