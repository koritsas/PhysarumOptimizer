package edu.koritsas.slimemold.shortestpath;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
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
 * Created by ilias on 6/11/2016.
 */
public abstract class PhysarumPolycephalumLangrarianCSP extends PhysarumPolycephalumSP {
    private double λ=0;
    private double step;
    private double λmax;
    private HashMap<Edge,Double> L;

    public PhysarumPolycephalumLangrarianCSP(Graph graph, Node sourceNode, Node sinkNode, double absoluteThreshold, double relativeThreshold, int numberOfIterations,double λmax,double step) {
        super(graph, sourceNode, sinkNode, absoluteThreshold, relativeThreshold, numberOfIterations);
        this.λmax=λmax;
        this.step=step;
    }


    @Override
    protected void initializeMaps(Graph graph) {
        super.initializeMaps(graph);
        Collection<Edge> edges =graph.getEdges();
        L=new HashMap<>(graph.getEdges().size());

        for (Edge e:edges){
            L.putIfAbsent(e,getEdgeCost(e));

        }
    }

    @Override
    protected double calculateCoefficient(Node n1, Node n2) {
        Edge e=n1.getEdge(n2);
        double coeff=0;
        if (e!=null){
            double D =conductivityMap.get(e);
            double w = L.get(e);
            coeff=-D/w;
        }

        return coeff;
    }

    @Override
    public void execute() {

        initializeMaps(graph);
        logger.info("Starting iterations...");
    while (λ<λmax) {
        // for (int i = 0; i <numberOfIterations ; i++) {

        List<Edge> edges = new ArrayList<>(graph.getEdges());
        edges.stream().forEach(edge -> L.put(edge, L.get(edge) + λ * getEdgeConstraintValue(edge)));
        while (!converged) {




            logger.info("Current iteration: " + iteration);
            List<Node> allNodes = getAllNodes(graph);


            List<Node> allButSinkNodes = getAllButSink(graph);


            RealVector constants = createConstantsVector();

            RealMatrix coefficients = createCoefficientsMatrix();


            // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();

            RealVector solution = solver.solve(constants);

            redefinePressures(allButSinkNodes, solution);

            redefineFlows(graph);

            redefineDiameters(graph);

            converged = checkConvergenceCriterion(previousFluxMap, currentFluxMap, iteration, m_absoluteThreshold, m_relativeThreshold);


            previousFluxMap = new HashMap<>(currentFluxMap);
            iteration++;
            if (iteration==numberOfIterations){
                break;
            }

        }
        if (iteration==numberOfIterations){
            break;
        }
        converged=false;
        boolean conv =pathViolatesConstraints(getGraph());
        if(conv==false){
            break;
        }

        double C= (double) getGraph().getEdges().stream().collect(Collectors.summingDouble(new ToDoubleFunction<Edge>() {
            @Override
            public double applyAsDouble(Edge value) {
                return getEdgeConstraintValue(value);
            }
        }));
       if (C<100){
           break;
       }
        λ=λ+step;
        }

    }
    public abstract boolean pathViolatesConstraints(Graph graph);

    public abstract double getEdgeConstraintValue(Edge edge);
}
