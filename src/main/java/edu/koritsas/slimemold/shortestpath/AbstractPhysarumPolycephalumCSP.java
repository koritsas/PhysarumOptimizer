package edu.koritsas.slimemold.shortestpath;

import edu.koritsas.slimemold.AbstractPhysarumPolycephalum;
import org.geotools.factory.Hints;
import org.geotools.graph.path.Path;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Created by ilias on 31/10/2016.
 */
public abstract class AbstractPhysarumPolycephalumCSP extends AbstractPhysarumPolycephalum {
    private double penaltyDivisor;
    private final  int  timesGrowingThreshold;
    private HashMap<Edge,Integer> growingTimes = new HashMap<>();
    private Path importantSegments= new Path();
    private boolean converged =false;

    public AbstractPhysarumPolycephalumCSP(Graph graph, Node sourceNode, Node sinkNode, double Io, double γ, int numberOfIterations,double penaltyDivisor,int timesGrowingThreshold) {
        super(graph, sourceNode, sinkNode, Io, γ, numberOfIterations);
        this.penaltyDivisor=penaltyDivisor;
        this.timesGrowingThreshold=timesGrowingThreshold;
    }

    @Override
    protected void initializeMaps(Graph graph) {
        super.initializeMaps(graph);
        growingTimes=new HashMap<Edge,Integer>(graph.getEdges().size());

    }

    @Override
    public void execute() {
        super.execute();
        judgePath(importantSegments);


    }
    private void judgePath(Path importantPath){
        if(importantPath.isValid()){

            if (calculatePathCost(importantPath)>calculateMaximumConstraintValue(graph)){

                List<Edge> edges = importantPath.getEdges();
                edges.stream().parallel().forEach(edge -> conductivityMap.put(edge,conductivityMap.get(edge)/penaltyDivisor));
                importantPath.clear();
            }else{
                converged=true;
            }






        }


    }

    private double calculatePathCost(Path path){
       double constraintValue= (double) path.stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
            @Override
            public double applyAsDouble(Edge edge) {

                return getEdgeConstraint(edge);
            }
        }));

      return constraintValue;

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

        if((newD-D)>1E-3){
            growingTimes.put(e,growingTimes.get(e)+1);
            if (growingTimes.get(e)>timesGrowingThreshold){

                importantSegments.addEdge(e);
                importantSegments.add(e.getNodeB());
                importantSegments.add(e.getNodeA());
            }
        }


        return newD;
    }

    public abstract double getEdgeConstraint(Edge e);

    public abstract double calculateMaximumConstraintValue(Graph graph);


}
