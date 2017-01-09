package edu.koritsas.slimemold;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.FastMath;
import org.geotools.graph.build.GraphBuilder;
import org.geotools.graph.build.basic.BasicDirectedGraphBuilder;
import org.geotools.graph.build.basic.BasicGraphBuilder;
import org.geotools.graph.structure.DirectedGraph;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.geotools.graph.structure.basic.BasicEdge;
import org.geotools.graph.structure.basic.BasicGraph;
import org.geotools.graph.structure.basic.BasicNode;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

import java.util.*;
import java.util.function.ToDoubleFunction;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Created by ilias on 5/10/2016.
 */
public abstract class AbstractPhysarumPolycephalum {
    protected Graph graph;

    public Node getSourceNode() {
        return sourceNode;
    }

    public void setSourceNode(Node sourceNode) {
        this.sourceNode = sourceNode;
    }

    public Node getSinkNode() {
        return sinkNode;
    }

    public void setSinkNode(Node sinkNode) {
        this.sinkNode = sinkNode;
    }

    protected Node sourceNode;
    protected Node sinkNode;
    protected List<Node> sinkNodesList;
    protected final double Io=100;

    protected  int numberOfIterations;
    protected HashMap<Edge,Double> currentFluxMap;
    protected HashMap<Edge,Double> previousFluxMap;
    protected HashMap<Edge,Double> conductivityMap;
    protected HashMap<Node,Double> pressureMap;
    protected boolean converged =false;

    protected final double m_absoluteThreshold;
    protected final double m_relativeThreshold;

    protected static final Logger logger = Logger.getLogger(AbstractPhysarumPolycephalum.class.getName());

    protected int iteration=0;

    protected HashMap<Edge, XYSeries> flowSeriesMap = new HashMap<Edge, XYSeries>();
    protected HashMap<Edge, XYSeries> conductivitySeriesMap =new HashMap<Edge, XYSeries>();
    protected HashMap<Node, XYSeries> pressureSeriesMap =new HashMap<Node, XYSeries>();

    public AbstractPhysarumPolycephalum(Graph graph, Node sourceNode, Node sinkNode,double absoluteThreshold,double relativeThreshold, int numberOfIterations){
        this.graph=graph;
        this.sourceNode=sourceNode;
        this.sinkNode=sinkNode;
        this.numberOfIterations=numberOfIterations;
        this.m_absoluteThreshold  = absoluteThreshold;
        this.m_relativeThreshold=relativeThreshold;

    }




    public HashMap<Edge, XYSeries> getFlowSeriesMap() {
        return flowSeriesMap;
    }

    public void setFlowSeriesMap(HashMap<Edge, XYSeries> flowSeriesMap) {
        this.flowSeriesMap = flowSeriesMap;
    }

    public HashMap<Edge, XYSeries> getConductivitySeriesMap() {
        return conductivitySeriesMap;
    }

    public void setConductivitySeriesMap(HashMap<Edge, XYSeries> conductivitySeriesMap) {
        this.conductivitySeriesMap = conductivitySeriesMap;
    }

    /**
     *
     * @return the graph
     */
    public Graph getGraph() {
        GraphBuilder builder=null;
       if (graph instanceof DirectedGraph){
           builder = new BasicDirectedGraphBuilder();
       }else if (graph instanceof BasicGraph){
           builder=new BasicGraphBuilder();
       }

       List<Edge> edges = (List<Edge>) graph.getEdges().stream().filter(o -> FastMath.abs(currentFluxMap.get(o))>=0.0000001).collect(Collectors.toList());
       for (Edge e:edges){
           builder.addEdge(e);
           if (!builder.getGraph().getNodes().contains(e.getNodeA())){
               builder.addNode(e.getNodeA());
           }
           if (!builder.getGraph().getNodes().contains(e.getNodeB())){
              builder.addNode(e.getNodeB());
           }

       }




        for (Edge e:edges){
          builder.addEdge(e);
        }



        return builder.getGraph();
    }



    /**
     * Executes the algorithm
     */
    public void execute() {

        double startingTime=System.currentTimeMillis();

        initializeMaps(graph);

        logger.info("Starting iterations...");

        while (converged==false){

            logger.info("Current iteration: "+iteration);

            List<Node> allButSinkNodes = getAllButSink(graph);

            RealVector constants = createConstantsVector();

            RealMatrix coefficients = createCoefficientsMatrix();



            // DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();

            RealVector solution = solver.solve(constants);

            redefinePressures(allButSinkNodes, solution);

            redefineFlows(graph);

            redefineDiameters(graph);

            converged =checkConvergenceCriterion(previousFluxMap,currentFluxMap,iteration,m_absoluteThreshold,m_relativeThreshold);

            previousFluxMap=new HashMap<>(currentFluxMap);

            iteration++;
        }
        double endingTime=System.currentTimeMillis();

        logger.info("Execution Time: "+ (startingTime-endingTime));
      //  eliminateEdges(graph);
    }

    protected boolean checkConvergenceCriterion(HashMap<Edge,Double> previousFluxMap,HashMap<Edge,Double> currentFluxMap, int currentIteration,double absoluteThreshold, double relativeThreshold){

        ConvergenceChecker<HashMap<Edge,Double>> checker = new ConvergenceChecker<HashMap<Edge, Double>>() {
            @Override
            public boolean converged(int iteration, HashMap<Edge, Double> previous, HashMap<Edge, Double> current) {

              boolean absolute=  previous.keySet().stream().allMatch(edge -> FastMath.abs(previous.get(edge)-current.get(edge))<absoluteThreshold);

               boolean relative =  previous.keySet().stream().allMatch(edge -> FastMath.max(FastMath.abs(previous.get(edge)), FastMath.abs(current.get(edge)))*relativeThreshold>=FastMath.abs(previous.get(edge)-current.get(edge)));

                boolean maxIterReached =currentIteration==numberOfIterations;

                current.values().stream().collect(Collectors.maxBy((o1, o2) -> 1));

                return absolute||relative||maxIterReached;
            }
        };




            return checker.converged(iteration,previousFluxMap,currentFluxMap);


    }
    /**
     *
     * @return the coefficients matrix A for the system Ax=B;
     */


    protected RealMatrix createCoefficientsMatrix() {

        List<Node> allNodes=getAllNodes(graph);
        List<Node> allButSink =getAllButSink(graph);

        double[][] matrix = new double[allNodes.size()][allButSink.size()];

        for (int i = 0; i < allNodes.size(); i++) {
            Node n1 = allNodes.get(i);
            for (int j = 0; j < allButSink.size(); j++) {
                Node n2 = allButSink.get(j);

                if(n1.equals(n2)){
                    matrix[i][j]=calculateSelfCoefficient(n1);
                }else{
                        matrix[i][j]=calculateCoefficient(n1,n2);
                }


            }
        }


        return new Array2DRowRealMatrix(matrix);
    }

    /**
     *
     * @return the constants vector B for the system Ax=B
     */

    protected RealVector createConstantsVector() {
        List<Node> allNodes =getAllNodes(graph);
        double[] constants = new double[allNodes.size()];
        for (int i = 0; i < constants.length; i++) {

            Node n1 = allNodes.get(i);
            if (n1.equals(sourceNode)) {

                    constants[i] = Io - calculateCoefficient(n1,sinkNode) * pressureMap.get(sinkNode);

            } else if (n1.equals(sinkNode)) {
                constants[i] = -Io - calculateSelfCoefficient(n1) *  pressureMap.get(sinkNode);
            } else {

                    constants[i] = -calculateCoefficient(n1,sinkNode) *  pressureMap.get(sinkNode);


            }
        }


        return new ArrayRealVector(constants);
    }
    /**
     * Initializes the XYSeries of the diagrams
     * @param graph
     */
    protected void initializeMaps(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        currentFluxMap =new HashMap<Edge, Double>(graph.getEdges().size());
        conductivityMap=new HashMap<Edge,Double>(graph.getEdges().size());
        for (Edge e:edges){
            currentFluxMap.put(e,0.0);
            conductivityMap.put(e,1.0);
        }

        Collection<Node> nodes =graph.getNodes();
        pressureMap=new HashMap<Node,Double>(nodes.size());
        for (Node n:nodes){
            pressureMap.put(n,0.0);

        }

        previousFluxMap=new HashMap<>(currentFluxMap);



    }



    /**
     * @param solution Redefines pressures on edges from the solution of the system
     */
    protected void redefinePressures(List<Node> allButSink, RealVector solution) {
        for (int i = 0; i < allButSink.size(); i++) {

            Node n = allButSink.get(i);
            //noinspection Since15
            pressureSeriesMap.putIfAbsent(n,new XYSeries("("+n.getID()+")"));
            pressureMap.put(n,solution.getEntry(i));
            pressureSeriesMap.get(n).add(iteration,solution.getEntry(i));
        }


    }

    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    @SuppressWarnings("Since15")
    protected void redefineFlows(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            flowSeriesMap.putIfAbsent(e,new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double Q =calculateTubeFlux(e);
            currentFluxMap.put(e,Q);
            flowSeriesMap.get(e).add(iteration,Q);
        }
       double sum=currentFluxMap.values().stream().collect(Collectors.summingDouble(new ToDoubleFunction<Double>() {
           @Override
           public double applyAsDouble(Double value) {
               return value;
           }
       }));


    }

    /**
     *
     * @return the cost of the solution
     */
    public double getSolutionCost() {

        Collection<Edge> edges = getGraph().getEdges();
        double cost = edges.stream().collect(Collectors.summingDouble(value -> getEdgeCost(value)));
        return  cost;
    }

    /**
     *
     * @param graph eliminates the edges of the graph that have conductivity near zero
     */
    protected void eliminateEdges(Graph graph){

        graph.getEdges().removeIf(o -> FastMath.abs(currentFluxMap.get(o))<0.0001);


    }


    /**
     *
     * @param graph
     * redefines the flows for the edges of the graph
     */
    protected void redefineDiameters(Graph graph){
        Collection<Edge> edges =graph.getEdges();
        for(Edge e:edges){
            //noinspection Since15
            conductivitySeriesMap.putIfAbsent(e, new XYSeries(e.getID()+"("+e.getNodeA().getID()+","+e.getNodeB().getID()+")"));
            double D =calculateTubeConductivity(e);
            conductivityMap.put(e,D);
            conductivitySeriesMap.get(e).add(iteration,D);
        }

    }

    /**
     *
     * @param e
     * @return the flux of the tube
     */
    protected double calculateTubeFlux(Edge e){

        double p1 =pressureMap.get(e.getNodeA());
        double p2 =pressureMap.get(e.getNodeB());

        double D =conductivityMap.get(e);
        double w = getEdgeCost(e);

        double Q=(D/w)*(p1-p2);


        return Q;
    }

    /**
     *
     * @param e
     * @return the new conductivity of the edge
     */
    public abstract double calculateTubeConductivity(Edge e);



    /**
     * show the pressure diagram
     */
    public void showPressureMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

/*
        Collection<Node> nodes =graph.getNodes();
        for(Node n:nodes){
            dataset.addSeries(pressureSeriesMap.get(n));
        }
    */
        Iterator it= pressureSeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart = ChartFactory.createXYLineChart("Time-Pressure", "Time", "H", dataset, PlotOrientation.VERTICAL,true,true,true);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");

        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }
    public void showConductivityMap(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );
        Iterator it= conductivitySeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }

        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Conductivity", "Time", "D", dataset, PlotOrientation.VERTICAL,true,false,false);
        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);


    }


    public void showFlowDiagram(){
        final XYSeriesCollection dataset = new XYSeriesCollection( );

        Iterator it= flowSeriesMap.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry) it.next();
            dataset.addSeries((XYSeries) pair.getValue());
            it.remove();
        }
        JFreeChart flowChart =ChartFactory.createXYLineChart("Time-Flow", "Time", "Q", dataset,PlotOrientation.VERTICAL,true,true,true);

        ChartPanel chartPanel = new ChartPanel( flowChart);
        chartPanel.setPreferredSize( new java.awt.Dimension( 560 , 367 ) );
        ApplicationFrame frame = new ApplicationFrame("dsafdsafda");
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    protected List<Node> getAllButSink(Graph graph) {




        return getAllNodes(graph).stream().filter(node -> node!=sinkNode).collect(Collectors.toList());
    }


    /**
     * @param graph
     * @return Returns a list with all the nodes in the graph
     */

    protected List<Node> getAllNodes(Graph graph) {


        List<Node> nodeList = new ArrayList<>(graph.getNodes());


        return nodeList;

    }

    /**
     * @param n
     * @return Returns the self-coefficient
     */

    public double calculateSelfCoefficient(Node n) {
        double selfC = 0;
        List<Edge> edges = n.getEdges();


       /* Iterator<Node> it =n.getRelated();

        while(it.hasNext()){
            selfC = selfC + FastMath.abs(calculateCoefficient(n,it.next()));
        }
*/
       for (Edge edge:edges){
           double D =conductivityMap.get(edge);
           double w = getEdgeCost(edge);
           selfC=selfC+FastMath.abs(-D/w);
       }

        return selfC;
    }

    protected double calculateCoefficient(Node n1, Node n2){

        //Edge e=n1.getEdge(n2);


       /* if (e!=null){
            double D =conductivityMap.get(e);
            double w = getEdgeCost(e);
            coeff=-D/w;
        }
*/
        List<Edge> edges = n1.getEdges(n2);
        double coeff=0;

        if (edges!=null){
            Iterator<Edge> edgeIterator =edges.iterator();
            while (edgeIterator.hasNext()){
                Edge edge=edgeIterator.next();
                double D =conductivityMap.get(edge);
                double w = getEdgeCost(edge);
                coeff=coeff-D/w;
            }

        }
        return coeff;
    }

    /**
     *
     * @param e
     * @return the weight of the edge e
     */
    public abstract double getEdgeCost(Edge e);







}
