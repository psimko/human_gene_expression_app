from dash import html, dcc

app_layout = html.Div([
    html.H1("Gene Expression - WHB-10Xv3 dataset", style={'textAlign': 'center'}),
    
    # Stores for precomputed figures (hidden storage)
    dcc.Store(id='store-fig-bin'),
    dcc.Store(id='store-fig-trin'),
    dcc.Store(id='store-genes'),
    # dcc.Store(id='store-expression-div'),
    # dcc.Store(id='store-expression-class'),
    dcc.Store(id='store-expression-neuronal-cells'),
    dcc.Store(id='store-expression-nonneuronal-cells'),
    # dcc.Store(id='store-expression-subclass'),
    # dcc.Store(id='store-expression-supertype'),

    html.Label("Enter genes (comma-separated):"),
    dcc.Input(id='gene-input', type='text', value='Aif1, Gfap', debounce=True, style={'width': '100%'}),

    html.Div([
    html.Button('Update Plot', id='update-button', n_clicks=0)
    ], id="button-container"),
    
    html.Div([  
        # First Sunburst Plot
        html.Div([
            dcc.Graph(id='sunburst1', style={'height': '600px'})  # Fixed height for consistency
        ], style={'display': 'inline-block', 'width': '48%', 'position': 'relative', 'verticalAlign': 'top'}),

        # Second Sunburst Plot with Toggle Overlayed
        html.Div([
            # Toggle Positioned Over the Graph (Absolute Positioning)
            html.Div([
                dcc.Checklist(
                    id='toggle-trinary',
                    options=[{'label': ' Show Trinarized Plot', 'value': 'trinary'}],
                    value=[],
                    style={'textAlign': 'left'}
                )
            ], style={'position': 'absolute', 'top': '0px', 'left': '0px', 'zIndex': '10'}),  

            dcc.Graph(id='sunburst2', style={'height': '600px'})  # Match first sunburst
        ], style={'display': 'inline-block', 'width': '48%', 'position': 'relative', 'verticalAlign': 'top'})

    ], style={'display': 'flex', 'justify-content': 'space-between', 'alignItems': 'center', 'marginTop': '20px', 'paddingBottom': '40px'}) #,

     #html.Div(id='gene-bar-plots', style={'marginTop': '20px'})
])
