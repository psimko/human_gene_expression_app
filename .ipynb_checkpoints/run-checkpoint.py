from app_context import app
from layout import app_layout
import callbacks  # Just to register callbacks
from local_config import PORT

app.layout = app_layout

# Callbacks automatically register from import if they use @app.callback or @callback

if __name__ == '__main__':
    app.run_server(host="0.0.0.0", port=PORT, debug=False)